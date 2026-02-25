import argparse
import json
import math
import os
import numpy as np
import pymol
from pymol import cmd

def vector_to_matrix(v1, v2):
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    dot = np.dot(v1, v2)
    if dot > 0.999999: return np.eye(4)
    if dot < -0.999999:
        M = np.eye(4); M[0,0], M[1,1], M[2,2] = -1, -1, 1
        return M
    cross = np.cross(v1, v2)
    s = np.linalg.norm(cross)
    vx = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - dot) / (s ** 2))
    M = np.eye(4); M[:3, :3] = R
    return M

def standardize_geometry(obj):
    """
    Standardizes receptor with Pocket at 0,0,0.
    Pore=Z, Receptor Center on -X (Pocket facing +X).
    """
    # 1. Pore to Z
    coords = cmd.get_coords(obj)
    evals, evecs = np.linalg.eigh(np.cov(coords.T))
    pore_vec = evecs[:, np.argmax(evals)]
    rot_m = vector_to_matrix(pore_vec, np.array([0.0, 0.0, 1.0]))
    cmd.transform_selection(obj, rot_m.flatten().tolist())
    
    # 2. Z-Flip (Ensure receptor bulk is mostly below the pocket-origin)
    z_center = np.mean(cmd.get_coords(obj), axis=0)[2]
    if z_center > 0:
        cmd.rotate("x", 180, obj, camera=0, origin=[0, 0, 0])

    # 3. X-Alignment (Rotate so Receptor COM is on -X axis)
    # This means the pocket (origin) is effectively facing +X
    com = np.mean(cmd.get_coords(obj), axis=0)
    angle = math.atan2(com[1], com[0])
    # Target angle for COM is 180 degrees (pi)
    rot_angle = math.pi - angle
    cmd.rotate("z", math.degrees(rot_angle), obj, camera=0, origin=[0, 0, 0])
    
    # Return the new COM offset
    final_com = np.mean(cmd.get_coords(obj), axis=0)
    return final_com.tolist()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True, help="Input receptor PDB")
    parser.add_argument("--json", required=True, help="Input metadata JSON")
    parser.add_argument("--out_pdb", required=True, help="Output aligned receptor PDB")
    parser.add_argument("--out_json", required=True, help="Output metadata JSON")
    parser.add_argument("--ref_pdb", help="Reference receptor PDB to align to")
    parser.add_argument("--ref_json", help="Reference metadata JSON")
    parser.add_argument("--master_receptor", help="Name of master reference receptor (e.g., '6X3U')")
    args = parser.parse_args()

    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()
    cmd.reinitialize()

    with open(args.json, 'r') as f:
        meta = json.load(f)

    # Load reference metadata if provided
    ref_meta = None
    if args.ref_json and os.path.exists(args.ref_json):
        with open(args.ref_json, 'r') as f:
            ref_meta = json.load(f)
    
    cmd.load(args.pdb, "target")
    resn = meta.get("target_ligand_resn")
    chain = meta.get("target_ligand_chain")
    lig_coords = cmd.get_coords(f"target and resn {resn} and chain {chain}") if resn else None

    # Initialize ligands dict in metadata
    ligands_dict = {}
    native_ligand_resn = None
    native_ligand_chain = None

    # Determine which alignment mode
    if lig_coords is not None and not args.ref_pdb:
        # Case A: Master Receptor - Independent alignment
        pocket_center = np.mean(lig_coords, axis=0)
        cmd.translate([-pocket_center[0], -pocket_center[1], -pocket_center[2]], "target")
        com_offset = standardize_geometry("target")
        meta["receptor_center_offset"] = com_offset

        # Track native ligand (stays at Z unless ref_ligands provided)
        native_ligand_resn = resn
        native_ligand_chain = chain

    elif args.ref_pdb:
        # Case B/C: Align to Reference (holo or apo)
        cmd.load(args.ref_pdb, "ref")

        # 1. Global Pre-Alignment
        cmd.super("target", "ref")

        # 2. Alpha/Gamma Radial Interface Alignment
        cmd.pseudoatom("pocket_center", pos=[0.0, 0.0, 0.0])
        ref_sel = "ref and byres (ref within 10 of pocket_center)"

        bzd_chains = [c for c, desc in meta.get("chains", {}).items() if "_bzd" in desc]
        if bzd_chains:
            bzd_sel = " or ".join([f"chain {c}" for c in bzd_chains])
            target_sel = f"target and ({bzd_sel}) and byres (target within 10 of pocket_center)"
        else:
            target_sel = "target and byres (target within 10 of pocket_center)"

        # Execute targeted superposition
        cmd.super(target_sel, ref_sel)

        # 3. Calculate COM offset
        com_offset = np.mean(cmd.get_coords("target"), axis=0).tolist()
        meta["receptor_center_offset"] = com_offset

        # Track native ligand if present (will be bumped)
        if lig_coords is not None:
            native_ligand_resn = resn
            native_ligand_chain = chain

    else:
        # Fallback: Just Z-align if no reference or ligand
        standardize_geometry("target")
        meta["receptor_center_offset"] = None

    # Copy ligands from reference structure if aligning to reference
    if args.ref_pdb and cmd.count_atoms("ref") > 0 and ref_meta:
        # Parse reference metadata to get ligand chains
        ref_ligands = ref_meta.get("ligands", {})

        if ref_ligands:
            # If target has native ligand at Z, bump it first
            if native_ligand_resn and native_ligand_chain == 'Z':
                cmd.alter(f"target and chain Z and resn {native_ligand_resn}", "chain='Y'")
                cmd.sort()

                # Record native ligand position
                native_coords = cmd.get_coords("target and chain Y")
                if native_coords is not None:
                    native_com = np.mean(native_coords, axis=0).tolist()
                    ligands_dict['Y'] = {
                        "resname": native_ligand_resn,
                        "source": "native",
                        "native": True,
                        "position": native_com
                    }

            # Copy each ligand from reference based on metadata
            for chain_id, lig_info in ref_ligands.items():
                if cmd.count_atoms(f"ref and chain {chain_id}") > 0:
                    # Copy the chain
                    cmd.create("lig_copy", f"ref and chain {chain_id}")
                    cmd.create("target", "target or lig_copy")
                    cmd.delete("lig_copy")

                    # Record ligand info
                    lig_coords = cmd.get_coords(f"target and chain {chain_id}")
                    if lig_coords is not None:
                        lig_com = np.mean(lig_coords, axis=0).tolist()
                        ligands_dict[chain_id] = {
                            "resname": lig_info["resname"],
                            "source": lig_info.get("source", args.master_receptor),
                            "native": False,
                            "position": lig_com
                        }

    # If this is master receptor (no ref_pdb) and has native ligand, record it
    if not args.ref_pdb and native_ligand_resn:
        lig_coords = cmd.get_coords(f"target and chain {native_ligand_chain}")
        if lig_coords is not None:
            lig_com = np.mean(lig_coords, axis=0).tolist()
            ligands_dict[native_ligand_chain] = {
                "resname": native_ligand_resn,
                "source": "native",
                "native": True,
                "position": lig_com
            }

    # Update metadata
    if ligands_dict:
        meta["ligands"] = ligands_dict
    if args.master_receptor:
        meta["alignment_reference"] = args.master_receptor

    # Clean up old keys
    if "target_ligand_resn" in meta:
        del meta["target_ligand_resn"]
    if "target_ligand_chain" in meta:
        del meta["target_ligand_chain"]
    if "target_ligand_center" in meta:
        del meta["target_ligand_center"]

    # Save Output PDB
    cmd.save(args.out_pdb, "target")

    # Save updated metadata
    with open(args.out_json, 'w') as f:
        json.dump(meta, f, indent=4)

    cmd.quit()

if __name__ == "__main__":
    main()