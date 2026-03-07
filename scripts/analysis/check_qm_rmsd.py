#!/usr/bin/env python3
import sys
import argparse
import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from Bio.PDB import PDBParser

def main():
    parser = argparse.ArgumentParser(description="Calculate RMSD of QM-prepared ligand to a crystal reference")
    parser.add_argument('--sdf', required=True, type=str, help="Path to generated ligand SDF")
    parser.add_argument('--pdb', required=True, type=str, help="Path to receptor PDB")
    parser.add_argument('--metadata', required=True, type=str, help="Path to receptor metadata JSON")
    args = parser.parse_args()

    mol = next(Chem.SDMolSupplier(args.sdf, removeHs=False))

    with open(args.metadata) as f:
        meta = json.load(f)

    parser_pdb = PDBParser(QUIET=True)
    structure = parser_pdb.get_structure('receptor', args.pdb)

    ligand_res = None
    for model in structure:
        for chain in model:
            if chain.id == meta['target_ligand_chain']:
                for res in chain:
                    if res.resname == meta['target_ligand_resn']:
                        ligand_res = res
                        break

    if not ligand_res:
        print(f"Error: Could not find ligand {meta['target_ligand_resn']} in chain {meta['target_ligand_chain']}")
        sys.exit(1)

    coords = [atom.coord for atom in ligand_res if atom.element.strip() != 'H']

    ref_mol = Chem.RemoveHs(mol)
    ref_conf = ref_mol.GetConformer()
    
    if len(coords) != ref_mol.GetNumAtoms():
        print(f"Error: Atom count mismatch! Generated={ref_mol.GetNumAtoms()}, Crystal={len(coords)}")
        sys.exit(1)

    for i in range(len(coords)):
        x, y, z = float(coords[i][0]), float(coords[i][1]), float(coords[i][2])
        ref_conf.SetAtomPosition(i, (x, y, z))

    gen_no_h = Chem.RemoveHs(mol)
    rmsd = rdMolAlign.GetBestRMS(gen_no_h, ref_mol, prbId=0, refId=0)
    print(f"RMSD to crystal: {rmsd:.3f} Å")

if __name__ == "__main__":
    main()
