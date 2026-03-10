#!/usr/bin/env python3
"""
Check positions of crystal ligand vs docked pose
Simple coordinate analysis to understand what's happening.
"""

import sys
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser
from rdkit import Chem


def get_crystal_ligand_coords(pdb_path: Path, chain: str, resn: str):
    """Extract crystal ligand coordinates from PDB."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', str(pdb_path))

    coords = []
    atoms = []

    for model in structure:
        for chain_obj in model:
            if chain_obj.id == chain:
                for residue in chain_obj:
                    if residue.resname == resn:
                        for atom in residue:
                            if atom.element != 'H':  # Skip hydrogens
                                coords.append(atom.coord)
                                atoms.append(atom.name)

    return np.array(coords), atoms


def get_docked_coords(sdf_path: Path):
    """Extract docked ligand coordinates from SDF."""
    mol = Chem.SDMolSupplier(str(sdf_path), removeHs=True)[0]

    if mol is None:
        return None, None

    conf = mol.GetConformer()
    coords = []
    atoms = []

    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        coords.append([pos.x, pos.y, pos.z])
        atoms.append(mol.GetAtomWithIdx(i).GetSymbol())

    return np.array(coords), atoms


def main():
    # Paths
    receptor_pdb = Path("data/receptors/6X3U_aligned.pdb")
    docked_sdf = Path("data/analysis/best_pose_conf1_pose14.sdf")

    print("="*60)
    print("LIGAND POSITION ANALYSIS")
    print("="*60)

    # Load crystal ligand
    crystal_coords, crystal_atoms = get_crystal_ligand_coords(receptor_pdb, "Z", "FYP")
    print(f"\nCrystal ligand (FYP, chain Z):")
    print(f"  Heavy atoms: {len(crystal_atoms)}")
    print(f"  Center: {crystal_coords.mean(axis=0)}")
    print(f"  Range X: [{crystal_coords[:,0].min():.2f}, {crystal_coords[:,0].max():.2f}]")
    print(f"  Range Y: [{crystal_coords[:,1].min():.2f}, {crystal_coords[:,1].max():.2f}]")
    print(f"  Range Z: [{crystal_coords[:,2].min():.2f}, {crystal_coords[:,2].max():.2f}]")

    # Load docked pose
    docked_coords, docked_atoms = get_docked_coords(docked_sdf)
    print(f"\nDocked pose (best RMSD):")
    print(f"  Heavy atoms: {len(docked_atoms)}")
    print(f"  Center: {docked_coords.mean(axis=0)}")
    print(f"  Range X: [{docked_coords[:,0].min():.2f}, {docked_coords[:,0].max():.2f}]")
    print(f"  Range Y: [{docked_coords[:,1].min():.2f}, {docked_coords[:,1].max():.2f}]")
    print(f"  Range Z: [{docked_coords[:,2].min():.2f}, {docked_coords[:,2].max():.2f}]")

    # Calculate separation
    crystal_center = crystal_coords.mean(axis=0)
    docked_center = docked_coords.mean(axis=0)
    separation = np.linalg.norm(crystal_center - docked_center)

    print(f"\n{'='*60}")
    print(f"Center-to-center distance: {separation:.3f} Å")
    print(f"{'='*60}")

    # If they match atom counts, calculate RMSD
    if len(crystal_coords) == len(docked_coords):
        # Simple RMSD (no alignment)
        diff = crystal_coords - docked_coords
        rmsd_noalign = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        print(f"RMSD (no alignment): {rmsd_noalign:.3f} Å")

        # Check if docking box was centered on crystal ligand
        print(f"\nDocking box center was [0, 0, 0]")
        print(f"Crystal ligand center: {crystal_center}")
        print(f"→ Box was {'CORRECTLY' if np.linalg.norm(crystal_center) < 5 else 'NOT'} centered on crystal ligand")
    else:
        print(f"\nWARNING: Atom count mismatch!")
        print(f"  Crystal: {len(crystal_coords)} atoms")
        print(f"  Docked: {len(docked_coords)} atoms")

    print(f"{'='*60}\n")


if __name__ == "__main__":
    main()
