import argparse
import json
import math
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
    if z_center > 0: cmd.rotate("x", 180, obj)
    
    # 3. X-Alignment (Rotate so Receptor COM is on -X axis)
    # This means the pocket (origin) is effectively facing +X
    com = np.mean(cmd.get_coords(obj), axis=0)
    angle = math.atan2(com[1], com[0])
    # Target angle for COM is 180 degrees (pi)
    rot_angle = math.pi - angle
    cmd.rotate("z", math.degrees(rot_angle), obj)
    
    # Return the new COM offset
    final_com = np.mean(cmd.get_coords(obj), axis=0)
    return final_com.tolist()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--out_pdb", required=True)
    parser.add_argument("--out_json", required=True)
    parser.add_argument("--ref_pdb", help="Standardized Flumazenil PDB (Pocket at 0,0,0)")
    args = parser.parse_args()

    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()
    cmd.reinitialize()
    
    with open(args.json, 'r') as f:
        meta = json.load(f)
    
    cmd.load(args.pdb, "target")
    resn = meta.get("target_ligand_resn")
    chain = meta.get("target_ligand_chain")
    lig_coords = cmd.get_coords(f"target and resn {resn} and chain {chain}") if resn else None

    if lig_coords is not None:
        # Case A: Holo - Move pocket to 0,0,0 and standardize
        pocket_center = np.mean(lig_coords, axis=0)
        cmd.translate([-pocket_center[0], -pocket_center[1], -pocket_center[2]], "target")
        com_offset = standardize_geometry("target")
        meta["receptor_center_offset"] = com_offset
    elif args.ref_pdb:
        # Case B: Apo - Snap to Reference
        cmd.load(args.ref_pdb, "ref")
        # Align target to reference. Since ref has pocket at 0,0,0, target now does too.
        cmd.super("target", "ref")
        com_offset = np.mean(cmd.get_coords("target"), axis=0).tolist()
        meta["receptor_center_offset"] = com_offset
    else:
        # Fallback: Just Z-align if no reference or ligand
        standardize_geometry("target")
        meta["receptor_center_offset"] = None

    # Save Output PDB
    cmd.save(args.out_pdb, "target")
    
    # Update Metadata
    if "target_ligand_center" in meta: del meta["target_ligand_center"]
    with open(args.out_json, 'w') as f:
        json.dump(meta, f, indent=4)
        
    # Full Receptor Visualization
    cmd.hide("all")
    cmd.show("cartoon", "target")
    cmd.view() 
    cmd.turn('y', 90) # Side view of pore
    cmd.png(args.out_pdb.replace(".pdb", ".png"), width=1200, height=1200, ray=1)
    
    cmd.quit()

if __name__ == "__main__":
    main()