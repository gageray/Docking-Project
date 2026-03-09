#!/usr/bin/env python3
"""
Convert CREST conformer ensemble (XYZ) to multi-structure SDF file.
Preserves all conformers from the CREST output.

Usage:
    python crest_ensemble_to_sdf.py -i crest_conformers.xyz -o output.sdf --smiles "SMILES_STRING"
"""

import argparse
import sys
import zipfile
import io
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.ML.Cluster import Butina

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


def xyz_to_sdf(xyz_file, output_sdf, smiles, prune=False, rmsd_thresh=1.0):
    """Convert CREST XYZ ensemble to multi-structure SDF.

    Args:
        xyz_file: Path to crest_conformers.xyz
        output_sdf: Path to output SDF file
        smiles: SMILES string for the molecule (defines connectivity)
    """
    # Parse XYZ file
    conformers = parse_xyz_ensemble(xyz_file)
    print(f"[*] Parsed {len(conformers)} conformers from {xyz_file}")

    # Create RDKit molecule from SMILES
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
    out_path = Path(output_sdf)
    if out_path.suffix != '.zip':
        out_path = out_path.with_suffix('.zip')

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

            # Write to string buffer
            sio = io.StringIO()
            writer = Chem.SDWriter(sio)
            writer.write(out_mol)
            writer.close()

            # Add to zip
            sdf_filename = f"{out_path.stem}_conf_{cid+1}.sdf"
            zf.writestr(sdf_filename, sio.getvalue())
            print(f"  ✓ Zipped {sdf_filename}: {energy:.4f} kcal/mol")

    print(f"\n[*] Export complete: {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Convert CREST XYZ ensemble to SDF")
    parser.add_argument('-i', '--input', required=True, help='Input crest_conformers.xyz file')
    parser.add_argument('-o', '--output', required=True, help='Output SDF file')
    parser.add_argument('-s', '--smiles', required=True, help='SMILES string (defines connectivity)')
    parser.add_argument('-p', '--prune', action='store_true', help='Enable Butina clustering to prune redundant conformers')
    parser.add_argument('-r', '--rmsd', type=float, default=1.0, help='RMSD threshold for pruning')

    args = parser.parse_args()

    xyz_to_sdf(args.input, args.output, args.smiles, args.prune, args.rmsd)


if __name__ == "__main__":
    main()
