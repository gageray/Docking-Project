import os
import json
import pymol
from pymol import cmd
import subprocess
import glob

def fix_pdb_elements(pdb_file):
    """
    RDKit/Meeko requires element symbols in columns 77-78.
    PDB2PQR often leaves these blank. This script restores them.
    """
    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    fixed_lines = []
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Atom name is in columns 12-16. 
            # Usually ' N  ', ' CA ', ' C  ', ' O  ', ' H  '
            atom_name = line[12:16].strip()
            
            # Extract the element: 
            # 1. If it starts with a digit (e.g., 1H), element is the second char
            # 2. Otherwise, it's the first char (N, C, O, S, P)
            element = atom_name[0] if not atom_name[0].isdigit() else atom_name[1]
            
            # Rebuild the line with the element symbol in columns 77-78
            # Ensure the line is long enough (80 chars)
            line = line.rstrip().ljust(76)
            line = line[:76] + element.rjust(2) + "\n"
        fixed_lines.append(line)

    with open(pdb_file, 'w') as f:
        f.writelines(fixed_lines)

def process_and_convert(pdb_path, meta_path, out_dir):
    """
    Load a prepped PDB, strip to BZD chains only (no heteroatoms/water),
    and convert to PDBQT.
    """
    if not os.path.exists(pdb_path) or not os.path.exists(meta_path):
        print(f"Skipping {pdb_path} or missing metadata at {meta_path}")
        return

    print(f"[*] Processing {pdb_path}...")
    with open(meta_path, 'r') as f:
        meta = json.load(f)

    # Clean PyMOL
    cmd.reinitialize()
    cmd.load(pdb_path, "receptor")

    # Keep only the BZD chains as defined in metadata
    bzd_chains = [ch for ch, name in meta.get("chains", {}).items() if "_bzd" in name]
    if not bzd_chains:
        print(f"[-] WARNING: No BZD chains found in metadata for {pdb_path}. Skipping.")
        return

    chain_sel = " or ".join([f"chain {c}" for c in bzd_chains])
    
    # Strip down to only the relevant BZD protein chains
    cmd.remove(f"not ({chain_sel})")
    
    # Strictly remove all lipids, ligands, and waters regardless of HETATM tags
    cmd.remove("not polymer.protein")

    # Name the output files
    base_name = os.path.basename(pdb_path).replace("_ligand", "").replace("_apo", "").replace("_prepped", "").split('.')[0]
    stripped_pdb = os.path.join(out_dir, f"{base_name}_stripped.pdb")
    final_pdbqt = os.path.join(out_dir, f"{base_name}_apo.pdbqt")

    # Save the stripped, unprotonated PDB temporarily
    raw_pdb = stripped_pdb.replace(".pdb", "_raw.pdb")
    cmd.save(raw_pdb, "receptor")
    print(f"[*] Saved raw stripped PDB to {raw_pdb}")

    # Add chemically accurate hydrogens at pH 7.4 using PDB2PQR
    # --keep-chain prevents it from stripping chain IDs (critical for GNINA tracking)
    protonated_pdb = stripped_pdb
    pdb2pqr_cmd = [
        "pdb2pqr",
        "--ff=PARSE",
        "--titration-state-method=propka",
        "--with-ph=7.4",
        "--keep-chain",
        "--nodebump",
        raw_pdb,
        protonated_pdb
    ]
    
    try:
        subprocess.run(pdb2pqr_cmd, check=True, capture_output=True)
        print(f"[+] Successfully protonated via PDB2PQR to {protonated_pdb}")
        
        # SANITIZE: Fix missing element columns for RDKit/Meeko
        fix_pdb_elements(protonated_pdb)
        print(f"[*] Sanitized element symbols in {protonated_pdb}")
        
    except subprocess.CalledProcessError as e:
        print(f"[-] PDB2PQR failed for {raw_pdb}:")
        if e.stderr:
            print(e.stderr.decode('utf-8'))
        return

    # Convert the accurately protonated PDB to PDBQT using Meeko
    # CRITICAL: Use --read_pdb (NOT -i) to bypass ProDy coordinate normalization
    # ProDy's "smart" behavior re-centers coordinates, breaking our alignment
    out_basename = os.path.splitext(final_pdbqt)[0]

    meeko_cmd = [
        "mk_prepare_receptor.py",
        "--read_pdb", protonated_pdb,  # Raw PDB reader - preserves coordinates
        "-o", out_basename,
        "-p",  # Write PDBQT output
        "--charge_model", "gasteiger"  # Gasteiger charges for GNINA scoring
    ]
    
    try:
        subprocess.run(meeko_cmd, check=True, capture_output=True)
        print(f"[+] Successfully converted via Meeko to {final_pdbqt}")
    except subprocess.CalledProcessError as e:
        print(f"[-] Meeko conversion failed for {protonated_pdb}:")
        if e.stderr:
            print(e.stderr.decode('utf-8'))
        
    # Clean up PyMOL for the next run
    cmd.delete("all")

import argparse

def main():
    parser = argparse.ArgumentParser(description="Strip and convert receptors to PDBQT format.")
    parser.add_argument("--pdb", type=str, required=True, help="Path to the prepped input PDB file.")
    parser.add_argument("--meta", type=str, required=True, help="Path to the corresponding metadata JSON file.")
    parser.add_argument("--out-dir", type=str, required=True, help="Directory to save the stripped PDB and output PDBQT.")
    args = parser.parse_args()

    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()
    
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir, exist_ok=True)
        
    process_and_convert(args.pdb, args.meta, args.out_dir)

    cmd.quit()
    print("[*] Conversion process complete.")

if __name__ == "__main__":
    main()
