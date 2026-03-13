#!/usr/bin/env python3
"""
Convert CREST conformer ensemble (XYZ) to PDBQT zip archive.
Preserves all conformers from CREST with aromatic atom types via Meeko.

Usage:
    python crest_ensemble_to_mol2.py -i crest_conformers.xyz -o output.zip --smiles "SMILES_STRING"
"""

import argparse
import sys
import zipfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.ML.Cluster import Butina

try:
    from meeko import MoleculePreparation
    MEEKO_AVAILABLE = True
except ImportError:
    MEEKO_AVAILABLE = False
    print("ERROR: Meeko not available. Install with: conda install -c conda-forge meeko")
    sys.exit(1)

def prune_redundant_conformers(ensemble_mol, rmsd_thresh):
    """
    Aligns conformers, calculates RMSD matrix, and uses Butina clustering 
    to drop redundant conformers. Returns list of conformer IDs to keep.
    """
    rdMolAlign.AlignMolConformers(ensemble_mol)
    rmsd_matrix = AllChem.GetConformerRMSMatrix(ensemble_mol, prealigned=True)
    
    # Butina clustering returns tuples of clustered IDs. Index 0 is the centroid.
    clusters = Butina.ClusterData(rmsd_matrix, ensemble_mol.GetNumConformers(), rmsd_thresh, isDistData=True)
    
    keep_ids = [cluster[0] for cluster in clusters]
    keep_ids.sort()
    return keep_ids

def parse_xyz_ensemble(xyz_file):
    """Parse CREST XYZ file containing multiple conformers.

    Returns:
        List of (energy, coords, elements) tuples
    """
    with open(xyz_file) as f:
        lines = f.readlines()

    conformers = []
    i = 0

    while i < len(lines):
        # First line: number of atoms
        natoms = int(lines[i].strip())

        # Second line: energy
        energy_line = lines[i+1].strip()
        try:
            energy = float(energy_line)
        except:
            energy = 0.0

        # Next natoms lines: element x y z
        coords = []
        elements = []
        for j in range(i+2, i+2+natoms):
            parts = lines[j].split()
            elements.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

        conformers.append((energy, coords, elements))
        i += natoms + 2

    return conformers


def xyz_to_pdbqt(xyz_file, output_zip, smiles, prune=False, rmsd_thresh=1.0):
    """Convert CREST XYZ ensemble to PDBQT zip archive via Meeko.

    Args:
        xyz_file: Path to crest_conformers.xyz
        output_zip: Path to output zip file (will contain PDBQT files)
        smiles: SMILES string (defines connectivity and aromaticity)
    """
    # Parse XYZ file
    conformers = parse_xyz_ensemble(xyz_file)
    print(f"[*] Parsed {len(conformers)} conformers from {xyz_file}")

    # Build template strictly from SMILES to map aromaticity and connectivity
    # Uses same RDKit logic as qm_ligand_prep.py, ensuring identical atom indexing
    mol_template = Chem.MolFromSmiles(smiles)
    if mol_template is None:
        print(f"ERROR: Invalid SMILES: {smiles}")
        sys.exit(1)
    mol_template = Chem.AddHs(mol_template)

    print(f"[*] Template molecule: {mol_template.GetNumAtoms()} atoms")

    # Verify atom count matches
    if len(conformers[0][2]) != mol_template.GetNumAtoms():
        print(f"ERROR: Atom count mismatch!")
        print(f"  XYZ file: {len(conformers[0][2])} atoms")
        print(f"  SMILES: {mol_template.GetNumAtoms()} atoms")
        sys.exit(1)

    # Build ensemble molecule to compute distance matrix
    ensemble_mol = Chem.Mol(mol_template)
    for idx, (energy, coords, elements) in enumerate(conformers):
        conf = Chem.Conformer(ensemble_mol.GetNumAtoms())
        for atom_idx, (x, y, z) in enumerate(coords):
            conf.SetAtomPosition(atom_idx, (x, y, z))
        ensemble_mol.AddConformer(conf, assignId=True)

    if prune:
        print(f"[*] Pruning conformers with RMSD < {rmsd_thresh} Å...")
        keep_ids = prune_redundant_conformers(ensemble_mol, rmsd_thresh)
        print(f"[*] Retaining {len(keep_ids)} of {len(conformers)} conformers")
    else:
        keep_ids = list(range(len(conformers)))

    # Ensure output has .zip extension
    out_path = Path(output_zip)
    if out_path.suffix != '.zip':
        out_path = out_path.with_suffix('.zip')

    # Initialize Meeko preparator
    # Default Meeko will add rotatable bonds, but GNINA can handle QM conformers
    preparator = MoleculePreparation()

    with zipfile.ZipFile(out_path, 'w', zipfile.ZIP_DEFLATED) as zf:
        for cid in keep_ids:
            energy, coords, elements = conformers[cid]

            out_mol = Chem.Mol(mol_template)
            out_mol.RemoveAllConformers()

            conf = Chem.Conformer(out_mol.GetNumAtoms())
            for atom_idx, (x, y, z) in enumerate(coords):
                conf.SetAtomPosition(atom_idx, (x, y, z))
            out_mol.AddConformer(conf)

            out_mol.SetProp("_Name", f"conformer_{cid+1}")
            out_mol.SetProp("CREST_Energy", f"{energy:.6f}")
            out_mol.SetProp("Conformer_ID", str(cid+1))

            # Pass RDKit Mol directly to Meeko (preserves aromatic flags from SMILES)
            try:
                preparator.prepare(out_mol)
                pdbqt_content = preparator.write_pdbqt_string()

                # Write PDBQT to zip
                pdbqt_filename = f"{out_path.stem}_conf_{cid+1}.pdbqt"
                zf.writestr(pdbqt_filename, pdbqt_content)
                print(f"  ✓ Zipped {pdbqt_filename}: {energy:.4f} kcal/mol")

            except Exception as e:
                print(f"  ✗ Meeko failed for conformer {cid+1}: {e}")

    print(f"\n[*] Export complete: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Convert CREST XYZ ensemble to PDBQT zip archive via Meeko")
    parser.add_argument('-i', '--input', required=True, help='Input crest_conformers.xyz file')
    parser.add_argument('-o', '--output', required=True, help='Output zip file (will contain PDBQT files)')
    parser.add_argument('-s', '--smiles', required=True, help='SMILES string (defines connectivity and aromaticity)')
    parser.add_argument('-p', '--prune', action='store_true', help='Enable Butina clustering to prune redundant conformers')
    parser.add_argument('-r', '--rmsd', type=float, default=1.0, help='RMSD threshold for pruning (default: 1.0 Å)')

    args = parser.parse_args()

    xyz_to_pdbqt(args.input, args.output, args.smiles, args.prune, args.rmsd)


if __name__ == "__main__":
    main()
