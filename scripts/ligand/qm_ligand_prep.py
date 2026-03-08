#!/usr/bin/env python3
"""
Single-Job QM Ligand Preparation Pipeline
Generates 3D coordinates and executes GFN2-xTB conformational searches.
Operates via CLI arguments or a single-job JSON config.

Uses QupKake QM-based pKa prediction to determine correct protonation state at target pH.
"""

import sys
import os
import json
import argparse
import subprocess
import shutil
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, Descriptors

# Import QupKake-based pKa prediction
try:
    from qupkake_pka import get_pka_and_state
    QUPKAKE_AVAILABLE = True
except ImportError as e:
    QUPKAKE_AVAILABLE = False
    print(f"WARNING: qupkake_pka module not found. Error: {e}")

    # Fallback: no-op function that returns input SMILES and empty pKa dict
    def get_pka_and_state(input_smiles: str, target_ph: float = 7.4, name: str = "ligand"):
        """Fallback when QupKake is unavailable - returns input SMILES unchanged."""
        print(f"[{name}] QupKake unavailable, using input SMILES as-is")
        return input_smiles, {}

def find_longest_aliphatic_chain(mol: Chem.Mol) -> list:
    """
    Finds the longest contiguous chain of non-ring sp3 carbons using DFS.
    Returns an ordered list of atom indices.
    """
    valid_atoms = set()
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 6 and 
            not atom.IsInRing() and 
            atom.GetHybridization() == Chem.HybridizationType.SP3):
            valid_atoms.add(atom.GetIdx())
            
    # Build adjacency list
    adj = {a: [] for a in valid_atoms}
    for bond in mol.GetBonds():
        b, e = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if b in valid_atoms and e in valid_atoms:
            adj[b].append(e)
            adj[e].append(b)
            
    # DFS to find longest path
    longest_path = []
    
    def dfs(current, path):
        nonlocal longest_path
        if len(path) > len(longest_path):
            longest_path = list(path)
        for neighbor in adj[current]:
            if neighbor not in path:
                dfs(neighbor, path + [neighbor])
                
    for node in valid_atoms:
        dfs(node, [node])
        
    return longest_path

def run_crest_job(name: str, mol: Chem.Mol, workdir: Path, outdir: Path, threads: int,
                  formal_charge: int = 0, energy_window: float = 6.0, cinp_path: Path = None):
    """Executes the CREST binary, parses the output, and saves the SDF."""
    seed_xyz = workdir / "seed.xyz"
    conf = mol.GetConformer()

    with open(seed_xyz, 'w') as f:
        f.write(f"{mol.GetNumAtoms()}\n")
        f.write(f"CREST seed\n")
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol():2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")

    cmd = ['crest', seed_xyz.name, '--gfn2', '--chrg', str(formal_charge),
           '--ewin', str(energy_window), '--quick', '-T', str(threads),
           '--alpb', 'water']
    if cinp_path and cinp_path.exists():
        cmd.extend(['--cinp', cinp_path.name])
        
    print(f"[{name}] Executing GFN2-xTB conformational search...")
    crest_env = os.environ.copy()
    crest_env['OMP_NUM_THREADS'] = '1'
    crest_env['MKL_NUM_THREADS'] = '1'
    crest_env['OPENBLAS_NUM_THREADS'] = '1'
    
    process = subprocess.Popen(
        cmd,
        cwd=workdir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=crest_env
    )
    
    for line in process.stdout:
        print(line, end='')
        
    process.wait()

    if process.returncode != 0:
        print(f"[{name}] ERROR: CREST execution failed with code {process.returncode}")
        sys.exit(1)
        
    best_xyz = workdir / "crest_best.xyz"
    if not best_xyz.exists():
         print(f"[{name}] ERROR: Output crest_best.xyz not found.")
         sys.exit(1)
         
    with open(best_xyz) as f:
        lines = f.readlines()
        
    if len(lines) >= mol.GetNumAtoms() + 2:
        new_conf = Chem.Conformer(mol.GetNumAtoms())
        for i in range(mol.GetNumAtoms()):
            parts = lines[i+2].split()
            x, y, z = map(float, parts[1:4])
            new_conf.SetAtomPosition(i, (x, y, z))
            
        mol.RemoveAllConformers()
        mol.AddConformer(new_conf)
        
        out_sdf = outdir / f"{name}_qm.sdf"
        writer = Chem.SDWriter(str(out_sdf))
        writer.write(mol)
        writer.close()
        print(f"[{name}] Successfully exported QM-optimized ligand to {out_sdf}")

def process_ligand(args):
    """Generates the 3D seed and routes to standard or constrained QM search."""
    outdir = Path(args.outdir)
    workdir = outdir / f"{args.name}_crest"
    workdir.mkdir(parents=True, exist_ok=True)

    # Step 1: Get correct protonation state via QupKake (if pH specified)
    target_ph = getattr(args, 'ph', None)
    if target_ph is not None:
        print(f"[{args.name}] Running QupKake pKa prediction at pH {target_ph}")
        corrected_smiles, _ = get_pka_and_state(args.smiles, target_ph, args.name)
    else:
        print(f"[{args.name}] No pH specified, using input SMILES as-is (skipping pKa prediction)")
        corrected_smiles = args.smiles

    # Step 2: Parse corrected SMILES
    mol = Chem.MolFromSmiles(corrected_smiles)
    if mol is None:
        print(f"[{args.name}] ERROR: Invalid SMILES: {corrected_smiles}")
        sys.exit(1)

    # Get formal charge for CREST
    formal_charge = Chem.GetFormalCharge(mol)
    print(f"[{args.name}] Formal charge: {formal_charge:+d}")
        
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol, params) != 0:
        print(f"[{args.name}] ERROR: ETKDG embedding failed.")
        sys.exit(1)

    cinp_path = workdir / "constrain.cinp"

    if args.mode == 'constrained':
        tail_carbons = find_longest_aliphatic_chain(mol)
        
        if len(tail_carbons) < 4:
            print(f"[{args.name}] WARNING: Aliphatic chain too short ({len(tail_carbons)} carbons) for constraints. Running standard.")
            AllChem.MMFFOptimizeMolecule(mol)
            run_crest_job(args.name, mol, workdir, outdir, args.threads, formal_charge)
            return

        print(f"[{args.name}] Identified {len(tail_carbons)}-carbon aliphatic tail. Applying straight-chain constraints.")
        
        # 1. Straighten the tail
        conf = mol.GetConformer()
        for i in range(len(tail_carbons) - 3):
            idx1, idx2, idx3, idx4 = tail_carbons[i:i+4]
            rdMolTransforms.SetDihedralDeg(conf, idx1, idx2, idx3, idx4, 180.0)
            
        # 2. Collect hydrogens attached to the tail
        tail_atoms = set(tail_carbons)
        for c_idx in tail_carbons:
            atom = mol.GetAtomWithIdx(c_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    tail_atoms.add(neighbor.GetIdx())

        # 3. Relax geometry with MMFF while holding tail straight
        props = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=0)
        for atom_idx in tail_atoms:
            ff.MMFFAddPositionConstraint(atom_idx, 0.0, 1.0e4)
        ff.Initialize()
        ff.Minimize(maxIts=2000)
        
        # 4. Generate CREST constraint file (1-based indices)
        crest_indices = [str(i + 1) for i in tail_atoms]
        cinp_content = f"$constrain\n  force constant=1.0\n  atoms: {','.join(crest_indices)}\n$end\n"
        with open(cinp_path, 'w') as f:
            f.write(cinp_content)

        run_crest_job(args.name, mol, workdir, outdir, args.threads, formal_charge, cinp_path=cinp_path)

    else:
        # Standard Mode
        if mol.GetNumConformers() == 0:
            print(f"[{args.name}] WARNING: ETKDG failed earlier. Attempting bare minimum embedding.")
            AllChem.EmbedMolecule(mol)

        if mol.GetNumConformers() == 0:
            print(f"[{args.name}] ERROR: Could not generate 3D coordinates.")
            sys.exit(1)

        AllChem.MMFFOptimizeMolecule(mol)
        run_crest_job(args.name, mol, workdir, outdir, args.threads, formal_charge)

def main():
    parser = argparse.ArgumentParser(description="Single-Job QM Ligand Preparation with QupKake pKa Prediction")
    parser.add_argument('-c', '--config', type=str, help='Path to single-job JSON config file')
    parser.add_argument('-n', '--name', type=str, help='Ligand name')
    parser.add_argument('-s', '--smiles', type=str, help='Ligand SMILES (neutral form recommended)')
    parser.add_argument('-o', '--outdir', type=str, default='data/ligands', help='Output directory')
    parser.add_argument('-t', '--threads', type=int, default=4, help='CPU threads (default: 4)')
    parser.add_argument('-m', '--mode', type=str, choices=['standard', 'constrained'], default='standard', help='Sampling mode')
    parser.add_argument('--ph', type=float, default=None, help='Target pH for protonation state (if not specified, skips pKa prediction)')
    parser.add_argument('--pka-only', action='store_true', help='Only run QupKake pKa prediction and save results, skip CREST')

    args = parser.parse_args()

    # If config file is provided, load it and let CLI args override it
    if args.config:
        with open(args.config, 'r') as f:
            config_data = json.load(f)
        
        # Override args with config data if the arg was not explicitly set in CLI
        for key, value in config_data.items():
            if getattr(args, key) is None or (key in ['outdir', 'threads', 'mode'] and f"--{key}" not in sys.argv and f"-{key[0]}" not in sys.argv):
                setattr(args, key, value)

    # Validate required arguments
    if not args.name or not args.smiles:
        print("ERROR: --name and --smiles are required (via CLI or config).")
        parser.print_help()
        sys.exit(1)

    # pKa-only mode: Just run prediction and save results
    if args.pka_only:
        if not QUPKAKE_AVAILABLE:
            print("ERROR: --pka-only requires qupkake_pka module.")
            print("Ensure qupkake_pka.py is in the same directory and QupKake is installed.")
            sys.exit(1)

        target_ph = getattr(args, 'ph', None)
        if target_ph is None:
            print("ERROR: --pka-only requires --ph to be specified.")
            sys.exit(1)
        corrected_smiles, pka_records = get_pka_and_state(args.smiles, target_ph, args.name)

        # Parse and get charge
        mol = Chem.MolFromSmiles(corrected_smiles)
        if mol is None:
            print(f"[{args.name}] ERROR: Invalid corrected SMILES: {corrected_smiles}")
            sys.exit(1)

        formal_charge = Chem.GetFormalCharge(mol)

        # Save to text file
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        output_file = outdir / f"{args.name}_pka_prediction.txt"

        with open(output_file, 'w') as f:
            f.write(f"Ligand: {args.name}\n")
            f.write(f"Target pH: {target_ph}\n")
            f.write(f"\n")
            f.write(f"Input SMILES:\n{args.smiles}\n")
            f.write(f"\n")
            f.write(f"Corrected SMILES (pH {target_ph}):\n{corrected_smiles}\n")
            f.write(f"\n")
            f.write(f"Formal Charge: {formal_charge:+d}\n")
            f.write(f"\n")

            # Write predicted pKa values
            f.write(f"Predicted pKa Values (QupKake):\n")
            if not pka_records:
                f.write("  No ionizable sites detected.\n")
            else:
                for atom_info, val in pka_records.items():
                    f.write(f"  {atom_info}: {val:.2f}\n")
            f.write(f"\n")

            f.write(f"Heavy Atoms: {mol.GetNumAtoms()}\n")
            f.write(f"Molecular Weight: {Descriptors.MolWt(mol):.2f}\n")

        print(f"\n[{args.name}] pKa prediction results saved to: {output_file}")
        sys.exit(0)

    # Normal mode: Check for CREST and run full pipeline
    if not shutil.which('crest'):
        print("ERROR: CREST binary not found in PATH.")
        sys.exit(1)

    process_ligand(args)

if __name__ == "__main__":
    main()