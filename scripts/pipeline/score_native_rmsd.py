#!/usr/bin/env python3
"""
Calculate RMSD between docked ligand poses and native crystal structure ligand.
Standard approach for validating docking accuracy when native pose is available.
Uses RDKit for robust ligand RMSD calculation.
"""
import argparse
import sys
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import math

def calculate_ligand_rmsd(native_pdb, docked_sdf):
    """
    Calculate RMSD for all poses in docked SDF vs native crystal ligand.
    Returns list of (pose_num, rmsd, n_atoms) tuples.
    """
    results = []

    # Load native crystal structure ligand
    print(f"Loading native ligand: {native_pdb}")
    if native_pdb.endswith('.pdb'):
        native_mol = Chem.MolFromPDBFile(native_pdb, removeHs=False, sanitize=True)
    elif native_pdb.endswith('.sdf'):
        native_supplier = Chem.SDMolSupplier(native_pdb, removeHs=False, sanitize=True)
        native_mol = native_supplier[0] if len(native_supplier) > 0 else None
    else:
        print(f"Error: Unsupported file format {native_pdb}")
        return results

    if native_mol is None:
        print(f"Error: Could not load native ligand from {native_pdb}")
        return results

    native_heavy = sum(1 for atom in native_mol.GetAtoms() if atom.GetAtomicNum() > 1)
    print(f"Native ligand: {native_heavy} heavy atoms\n")

    # Load all docked poses from SDF
    print(f"Loading docked poses: {docked_sdf}")
    supplier = Chem.SDMolSupplier(docked_sdf, removeHs=False, sanitize=True)

    if supplier is None:
        print(f"Error: Could not load SDF file {docked_sdf}")
        return results

    num_poses = len(supplier)
    print(f"Docked poses: {num_poses} conformers\n")

    # Calculate RMSD for each pose
    for i, docked_mol in enumerate(supplier):
        pose_num = i + 1

        if docked_mol is None:
            print(f"Pose {pose_num}: Error - Could not read molecule")
            results.append((pose_num, float('nan'), 0))
            continue

        docked_heavy = sum(1 for atom in docked_mol.GetAtoms() if atom.GetAtomicNum() > 1)

        try:
            # Try best alignment first (handles symmetry if molecules match)
            rmsd = rdMolAlign.GetBestRMS(docked_mol, native_mol)
            n_atoms = min(docked_heavy, native_heavy)

            results.append((pose_num, rmsd, n_atoms))
            status = "EXCELLENT" if rmsd < 1.0 else ("GOOD" if rmsd < 2.0 else "POOR")
            print(f"Pose {pose_num}: RMSD = {rmsd:.3f} Å ({n_atoms} heavy atoms) [{status}]")

        except Exception as e:
            # If GetBestRMS fails (no substructure match), try CalcRMS (coordinate-only)
            try:
                rmsd = rdMolAlign.CalcRMS(docked_mol, native_mol)
                n_atoms = min(docked_heavy, native_heavy)

                results.append((pose_num, rmsd, n_atoms))
                status = "EXCELLENT" if rmsd < 1.0 else ("GOOD" if rmsd < 2.0 else "POOR")
                print(f"Pose {pose_num}: RMSD = {rmsd:.3f} Å ({n_atoms} heavy atoms, no alignment) [{status}]")
            except Exception as e2:
                print(f"Pose {pose_num}: Error - {e2}")
                results.append((pose_num, float('nan'), 0))

    return results

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate ligand RMSD between docked poses and native crystal structure"
    )
    parser.add_argument("--native", required=True, help="Native crystal ligand (PDB)")
    parser.add_argument("--docked", required=True, help="Docked poses (SDF)")
    parser.add_argument("--out", help="Output file (optional)")
    args = parser.parse_args()

    results = calculate_ligand_rmsd(args.native, args.docked)

    if not results:
        print("\nNo results obtained")
        sys.exit(1)

    # Find best pose
    valid_results = [(p, r, n) for p, r, n in results if not math.isnan(r)]
    if valid_results:
        best_pose, best_rmsd, best_atoms = min(valid_results, key=lambda x: x[1])
        print(f"\n{'='*50}")
        print(f"BEST POSE: {best_pose} with RMSD = {best_rmsd:.3f} Å")
        print(f"{'='*50}")

        # Write to output file if specified
        if args.out:
            with open(args.out, 'w') as f:
                f.write(f"Native: {args.native}\n")
                f.write(f"Docked: {args.docked}\n")
                f.write(f"\nPose\tRMSD (Å)\tAtoms\tQuality\n")
                for pose, rmsd, n_atoms in results:
                    if math.isnan(rmsd):
                        quality = "ERROR"
                    elif rmsd < 1.0:
                        quality = "EXCELLENT"
                    elif rmsd < 2.0:
                        quality = "GOOD"
                    else:
                        quality = "POOR"
                    f.write(f"{pose}\t{rmsd:.3f}\t{n_atoms}\t{quality}\n")
                f.write(f"\nBest: Pose {best_pose}, RMSD = {best_rmsd:.3f} Å\n")
            print(f"\nResults saved to {args.out}")
    else:
        print("\nAll poses failed to calculate RMSD")
        sys.exit(1)
