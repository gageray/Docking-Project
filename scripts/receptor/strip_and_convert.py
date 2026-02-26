import os
import json
import pymol
from pymol import cmd
import subprocess
import glob

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
    
    # Strip down to only the relevant BZD protein chains:
    # Remove everything that is not in the desired chains
    cmd.remove(f"not ({chain_sel})")
    
    # Remove heteroatoms (like the old ligand, crystallographic waters, etc.)
    # GNINA receptor should purely be the protein.
    cmd.remove("hetatm")
    cmd.remove("resn HOH")

    # Name the output files
    base_name = os.path.basename(pdb_path).replace("_ligand", "").replace("_apo", "").replace("_prepped", "").split('.')[0]
    stripped_pdb = os.path.join(out_dir, f"{base_name}_stripped.pdb")
    final_pdbqt = os.path.join(out_dir, f"{base_name}_apo.pdbqt")

    # Save the stripped PDB temporarily
    cmd.save(stripped_pdb, "receptor")
    print(f"[+] Saved temporary stripped PDB to {stripped_pdb}")

    # Convert to PDBQT using OpenBabel
    # -xr = strict atomic parsing? -p 7.4 adds protonation at phys pH.
    obabel_cmd = [
        "obabel", 
        "-i", "pdb", stripped_pdb, 
        "-o", "pdbqt", 
        "-O", final_pdbqt, 
        "-p", "7.4", 
        "--partialcharge", "gasteiger"
    ]
    
    try:
        subprocess.run(obabel_cmd, check=True, capture_output=True)
        print(f"[+] Successfully converted to {final_pdbqt}")
    except subprocess.CalledProcessError as e:
        print(f"[-] OpenBabel conversion failed for {stripped_pdb}:")
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
