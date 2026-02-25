import argparse
import sys
import numpy as np
import os
import json
import pymol
from pymol import cmd

def calc_rmsd(ref_file, target_file):
    cmd.reinitialize()
    cmd.load(ref_file, "ref")
    cmd.load(target_file, "target")
    
    cmd.pseudoatom("origin_pt", pos=[0,0,0])
    
    cmd.select("ref_pocket", "ref and name CA and byres (ref within 10.0 of origin_pt)")
    cmd.select("target_pocket", "target and name CA and byres (target within 10.0 of origin_pt)")
    
    try:
        align_str = cmd.super("target_pocket", "ref_pocket", cycles=0, transform=0)
        rmsd = align_str[0]
        n_atoms = align_str[1]
    except Exception as e:
        print(f"Error during alignment: {e}")
        rmsd = float('nan')
        n_atoms = 0
        
    return rmsd, n_atoms

def get_min_dist(target_file, lig_chain, res_chain, resi):
    cmd.reinitialize()
    cmd.load(target_file, "target")
    
    lig_sel = f"target and chain {lig_chain} and not hydro"
    res_sel = f"target and chain {res_chain} and resi {resi} and name CA"
    
    if cmd.count_atoms(lig_sel) == 0 or cmd.count_atoms(res_sel) == 0:
        return None
        
    dist = cmd.distance("tmp_dist", lig_sel, res_sel)
    cmd.delete("tmp_dist")
    return dist

if __name__ == "__main__":
    pymol.pymol_argv = ['pymol', '-qc']
    pymol.finish_launching()

    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True, help="Reference PDB file (6X3U_aligned.pdb)")
    parser.add_argument("--targets", nargs="+", required=True)
    parser.add_argument("--out", required=True, help="Output file for validation results")
    args = parser.parse_args()

    # Load reference metadata to get baseline distances
    ref_json = args.ref.replace("_aligned.pdb", "_aligned.json")
    if not os.path.exists(ref_json):
        print(f"Error: Reference metadata not found: {ref_json}")
        sys.exit(1)

    with open(ref_json, 'r') as f:
        ref_meta = json.load(f)

    # Get reference key residues
    ref_key_res = ref_meta.get("key_residues", {})
    ref_his_resi = ref_key_res.get("bzd_his", {}).get("resi", "102")  # Default to His102 for α1
    ref_phe_resi = ref_key_res.get("bzd_phe", {}).get("resi", "77")

    # Calculate baseline distances from reference FYP (chain Z)
    ref_dist_z_his = get_min_dist(args.ref, "Z", "D", ref_his_resi)
    ref_dist_z_phe = get_min_dist(args.ref, "Z", "E", ref_phe_resi)

    if ref_dist_z_his is None or ref_dist_z_phe is None:
        print(f"Error: Cannot calculate reference distances from {args.ref}")
        sys.exit(1)

    # Store output to file
    with open(args.out, "w") as f:
        f.write(f"Reference: {os.path.basename(args.ref)} (FYP baseline at chain Z)\n")
        f.write(f"Baseline distances: FYP-His{ref_his_resi}={ref_dist_z_his:.3f}Å, FYP-Phe{ref_phe_resi}={ref_dist_z_phe:.3f}Å\n")
        f.write("-" * 140 + "\n")
        f.write(f"{'Target':<20} | {'Pocket RMSD (CA)':<16} | {'CA Match':<10} | {'Native Lig':<12} | {'Δ FYP-His':<12} | {'Δ DZP-His':<12} | {'Δ FYP-Phe':<12} | {'Δ DZP-Phe':<12}\n")
        f.write("-" * 140 + "\n")

        for target in args.targets:
            # Load target metadata
            target_json = target.replace("_aligned.pdb", "_aligned.json")
            if not os.path.exists(target_json):
                print(f"Warning: Metadata not found for {target}, skipping")
                continue

            with open(target_json, 'r') as jf:
                target_meta = json.load(jf)

            # Find native ligand (if exists)
            native_resn = "APO"
            ligands = target_meta.get("ligands", {})
            for chain, lig_info in ligands.items():
                if lig_info.get("native", False):
                    native_resn = lig_info.get("resname", "???")
                    break

            # Calculate pocket RMSD
            rmsd, n_atoms = calc_rmsd(args.ref, target)

            # Use REFERENCE residues for all comparisons (to keep baseline consistent)
            # Distance from all ligands to key residues
            dist_z_his = get_min_dist(target, "Z", "D", ref_his_resi) if "Z" in ligands else None
            dist_z_phe = get_min_dist(target, "Z", "E", ref_phe_resi) if "Z" in ligands else None

            dist_y_his = get_min_dist(target, "Y", "D", ref_his_resi) if "Y" in ligands else None
            dist_y_phe = get_min_dist(target, "Y", "E", ref_phe_resi) if "Y" in ligands else None

            # Calculate deltas from REFERENCE FYP baseline (universal across all receptors)
            dz_his_str = f"{dist_z_his - ref_dist_z_his:+.3f}" if dist_z_his is not None else "N/A"
            dy_his_str = f"{dist_y_his - ref_dist_z_his:+.3f}" if dist_y_his is not None else "N/A"
            dz_phe_str = f"{dist_z_phe - ref_dist_z_phe:+.3f}" if dist_z_phe is not None else "N/A"
            dy_phe_str = f"{dist_y_phe - ref_dist_z_phe:+.3f}" if dist_y_phe is not None else "N/A"

            basename = os.path.basename(target)
            f.write(f"{basename:<20} | {rmsd:<16.3f} | {int(n_atoms):<10} | {native_resn:<12} | {dz_his_str:<12} | {dy_his_str:<12} | {dz_phe_str:<12} | {dy_phe_str:<12}\n")

    # Print to console as well
    with open(args.out, "r") as f:
        print(f.read())
    print(f"\nValidation metrics saved to {args.out}")
