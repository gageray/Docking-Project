import argparse
import sys
import numpy as np
import os
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
    parser.add_argument("--ref", required=True)
    parser.add_argument("--targets", nargs="+", required=True)
    parser.add_argument("--out", required=True, help="Output file for validation results")
    args = parser.parse_args()
    
    # First, calculate baseline distances from the reference structure (args.ref)
    ref_dist_z_his = get_min_dist(args.ref, "Z", "D", "101")
    if ref_dist_z_his is None: ref_dist_z_his = get_min_dist(args.ref, "Z", "D", "102")
    
    ref_dist_y_his = get_min_dist(args.ref, "Y", "D", "101")
    if ref_dist_y_his is None: ref_dist_y_his = get_min_dist(args.ref, "Y", "D", "102")

    ref_dist_z_phe = get_min_dist(args.ref, "Z", "E", "77")
    ref_dist_y_phe = get_min_dist(args.ref, "Y", "E", "77")
    
    # Store output to file
    with open(args.out, "w") as f:
        f.write(f"Reference: {os.path.basename(args.ref)}\n")
        f.write("-" * 115 + "\n")
        f.write(f"{'Target':<20} | {'Pocket RMSD (CA)':<16} | {'CA Match':<10} | {'Δ FYP(Z)-His102':<18} | {'Δ DZP(Y)-His102':<18}\n")
        f.write("-" * 115 + "\n")
        
        for target in args.targets:
            rmsd, n_atoms = calc_rmsd(args.ref, target)
            
            # Distance from Ligand Z (FYP) to alpha chain His101/102
            dist_z_his = get_min_dist(target, "Z", "D", "101")
            if dist_z_his is None: dist_z_his = get_min_dist(target, "Z", "D", "102") 
            
            # Distance from Ligand Y (DZP) to alpha chain His101/102
            dist_y_his = get_min_dist(target, "Y", "D", "101")
            if dist_y_his is None: dist_y_his = get_min_dist(target, "Y", "D", "102")
                
            dz_str = f"{dist_z_his - ref_dist_z_his:+.3f}" if (dist_z_his is not None and ref_dist_z_his is not None) else "N/A"
            # If DZP(Y) is present in target, compare its distance against reference FYP(Z) since 6X3U does not have a DZP validation ligand.
            dy_str = f"{dist_y_his - ref_dist_z_his:+.3f}" if (dist_y_his is not None and ref_dist_z_his is not None) else "N/A"
            
            basename = os.path.basename(target)
            f.write(f"{basename:<20} | {rmsd:<16.3f} | {int(n_atoms):<10} | {dz_str:<18} | {dy_str:<18}\n")
            
        f.write("\n")
        f.write("-" * 115 + "\n")
        f.write(f"{'Target':<20} | {'Pocket RMSD (CA)':<16} | {'CA Match':<10} | {'Δ FYP(Z)-Phe77':<18} | {'Δ DZP(Y)-Phe77':<18}\n")
        f.write("-" * 115 + "\n")
        
        for target in args.targets:
            rmsd, n_atoms = calc_rmsd(args.ref, target)
            
            # Distance from Ligand Z (FYP) to gamma chain Phe77 (chain E resi 77)
            dist_z_phe = get_min_dist(target, "Z", "E", "77")
            
            # Distance from Ligand Y (DZP) to gamma chain Phe77 
            dist_y_phe = get_min_dist(target, "Y", "E", "77")
                
            dz_str = f"{dist_z_phe - ref_dist_z_phe:+.3f}" if (dist_z_phe is not None and ref_dist_z_phe is not None) else "N/A"
            dy_str = f"{dist_y_phe - ref_dist_z_phe:+.3f}" if (dist_y_phe is not None and ref_dist_z_phe is not None) else "N/A"
            
            basename = os.path.basename(target)
            f.write(f"{basename:<20} | {rmsd:<16.3f} | {int(n_atoms):<10} | {dz_str:<18} | {dy_str:<18}\n")
            
    # Print to console as well
    with open(args.out, "r") as f:
        print(f.read())
    print(f"\nValidation metrics saved to {args.out}")
