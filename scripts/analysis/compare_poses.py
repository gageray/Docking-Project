#!/usr/bin/env python3
"""
Pose Comparison and Validation
Compares docked poses against reference crystal structure ligand using RMSD and interaction analysis.
"""

import sys
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from Bio.PDB import PDBParser, Selection


class PoseComparator:
    """Compares docked poses to reference crystal structure."""

    def __init__(self, receptor_pdb: Path, metadata_path: Path, docked_sdf: Path):
        """
        Initialize comparator.

        Args:
            receptor_pdb: Path to aligned receptor PDB with reference ligand
            metadata_path: Path to receptor metadata JSON
            docked_sdf: Path to GNINA output SDF with docked poses
        """
        self.receptor_pdb = Path(receptor_pdb)
        self.metadata_path = Path(metadata_path)
        self.docked_sdf = Path(docked_sdf)

        # Load metadata
        with open(self.metadata_path) as f:
            self.metadata = json.load(f)

        self.ref_ligand_resn = self.metadata['target_ligand_resn']
        self.ref_ligand_chain = self.metadata['target_ligand_chain']

        print(f"Target ligand: {self.ref_ligand_resn} chain {self.ref_ligand_chain}")

        # Extract reference ligand
        self.ref_mol = self._extract_reference_ligand()

        # Load docked poses
        self.docked_poses = self._load_docked_poses()

        print(f"Loaded {len(self.docked_poses)} docked poses")

    def _extract_reference_ligand(self) -> Chem.Mol:
        """Extract reference ligand from receptor PDB and convert to RDKit mol."""
        print(f"Extracting reference ligand from {self.receptor_pdb}")

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('receptor', str(self.receptor_pdb))

        # Find the ligand residue
        ligand_residue = None
        for model in structure:
            for chain in model:
                if chain.id == self.ref_ligand_chain:
                    for residue in chain:
                        if residue.resname == self.ref_ligand_resn:
                            ligand_residue = residue
                            break
                if ligand_residue:
                    break
            if ligand_residue:
                break

        if not ligand_residue:
            raise ValueError(f"Could not find ligand {self.ref_ligand_resn} in chain {self.ref_ligand_chain}")

        # Extract heavy atoms and coordinates
        atoms = []
        coords = []
        for atom in ligand_residue:
            # Skip hydrogens
            if atom.element.strip() != 'H':
                atoms.append({
                    'element': atom.element.strip(),
                    'name': atom.name,
                    'coord': atom.coord
                })
                coords.append(atom.coord)

        print(f"Extracted {len(atoms)} heavy atoms from reference ligand")

        # Create RDKit mol from coordinates
        # We'll use a simple approach: create mol from SDF in docked poses as template
        # and update coordinates to match PDB
        mol = Chem.RWMol()

        # For now, store as dict for RMSD calculation
        # We'll match to first docked pose to get proper RDKit mol
        self.ref_coords = np.array(coords)
        self.ref_atoms = atoms

        return None  # Will be set after loading first docked pose

    def _load_docked_poses(self) -> List[Tuple[Chem.Mol, Dict]]:
        """Load all docked poses from SDF file."""
        print(f"Loading docked poses from {self.docked_sdf}")

        poses = []
        supplier = Chem.SDMolSupplier(str(self.docked_sdf), removeHs=True, sanitize=True)

        for idx, mol in enumerate(supplier):
            if mol is None:
                print(f"WARNING: Could not parse molecule {idx}")
                continue

            # Extract GNINA properties
            props = {
                'pose_idx': idx,
                'gnina_score': float(mol.GetProp('minimizedAffinity')) if mol.HasProp('minimizedAffinity') else None,
                'cnn_score': float(mol.GetProp('CNNscore')) if mol.HasProp('CNNscore') else None,
                'cnn_affinity': float(mol.GetProp('CNNaffinity')) if mol.HasProp('CNNaffinity') else None,
            }

            poses.append((mol, props))

        return poses

    def calculate_rmsd(self, probe_mol: Chem.Mol, ref_mol: Optional[Chem.Mol] = None) -> float:
        """
        Calculate heavy atom RMSD between probe and reference.
        Uses symmetry-corrected alignment if possible.

        Args:
            probe_mol: Docked pose molecule
            ref_mol: Reference molecule (if None, uses self.ref_coords)

        Returns:
            RMSD in Angstroms
        """
        if ref_mol is not None:
            # Use RDKit's built-in RMSD with symmetry correction
            try:
                rmsd = rdMolAlign.GetBestRMS(probe_mol, ref_mol)
                return rmsd
            except Exception as e:
                print(f"WARNING: Symmetry-corrected RMSD failed: {e}, using simple RMSD")

        # Fallback: simple coordinate RMSD
        probe_conf = probe_mol.GetConformer()
        probe_coords = np.array([probe_conf.GetAtomPosition(i) for i in range(probe_mol.GetNumAtoms())])

        if probe_coords.shape[0] != self.ref_coords.shape[0]:
            print(f"WARNING: Atom count mismatch: probe={probe_coords.shape[0]}, ref={self.ref_coords.shape[0]}")
            return float('inf')

        # Calculate RMSD without alignment (assuming already aligned to same pocket)
        diff = probe_coords - self.ref_coords
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))

        return rmsd

    def compare_all_poses(self) -> List[Dict]:
        """
        Calculate RMSD for all docked poses.

        Returns:
            List of result dicts sorted by RMSD
        """
        print("Calculating RMSD for all poses...")

        results = []

        # Use first pose as reference mol template if needed
        if self.ref_mol is None and len(self.docked_poses) > 0:
            first_mol, _ = self.docked_poses[0]
            # Check atom counts match
            if first_mol.GetNumAtoms() == len(self.ref_coords):
                print("Reference structure matches docked pose atom count - using for template")
                # Create reference mol with correct coords
                self.ref_mol = Chem.Mol(first_mol)
                conf = self.ref_mol.GetConformer()
                for i in range(len(self.ref_coords)):
                    x, y, z = float(self.ref_coords[i][0]), float(self.ref_coords[i][1]), float(self.ref_coords[i][2])
                    conf.SetAtomPosition(i, (x, y, z))

        for mol, props in self.docked_poses:
            try:
                rmsd = self.calculate_rmsd(mol, self.ref_mol)

                result = {
                    'pose_idx': props['pose_idx'],
                    'rmsd': rmsd,
                    'gnina_score': props['gnina_score'],
                    'cnn_score': props['cnn_score'],
                    'cnn_affinity': props['cnn_affinity'],
                    'pass_threshold': rmsd < 2.0  # Standard threshold
                }

                results.append(result)

                print(f"Pose {props['pose_idx']}: RMSD={rmsd:.3f} Å, GNINA={props['gnina_score']:.3f}")

            except Exception as e:
                print(f"ERROR: Failed to calculate RMSD for pose {props['pose_idx']}: {e}")
                continue

        # Sort by RMSD
        results.sort(key=lambda x: x['rmsd'])

        return results

    def generate_report(self, results: List[Dict], output_path: Path):
        """Generate comparison report."""
        output_path = Path(output_path)

        # JSON report
        json_path = output_path.with_suffix('.json')
        with open(json_path, 'w') as f:
            json.dump({
                'receptor': str(self.receptor_pdb),
                'reference_ligand': f"{self.ref_ligand_resn}:{self.ref_ligand_chain}",
                'docked_sdf': str(self.docked_sdf),
                'total_poses': len(results),
                'passing_poses': sum(1 for r in results if r['pass_threshold']),
                'best_rmsd': results[0]['rmsd'] if results else None,
                'results': results
            }, f, indent=2)

        print(f"Saved JSON report to {json_path}")

        # CSV report
        csv_path = output_path.with_suffix('.csv')
        with open(csv_path, 'w') as f:
            f.write("pose_idx,rmsd,gnina_score,cnn_score,cnn_affinity,pass_threshold\n")
            for r in results:
                f.write(f"{r['pose_idx']},{r['rmsd']:.4f},{r['gnina_score']:.4f},"
                       f"{r['cnn_score']:.4f},{r['cnn_affinity']:.4f},{r['pass_threshold']}\n")

        print(f"Saved CSV report to {csv_path}")

        # Summary
        if results:
            print("\n" + "="*60)
            print("POSE COMPARISON SUMMARY")
            print("="*60)
            print(f"Total poses: {len(results)}")
            print(f"Passing poses (RMSD < 2.0 Å): {sum(1 for r in results if r['pass_threshold'])}")
            print(f"\nTop 5 poses by RMSD:")
            for i, r in enumerate(results[:5], 1):
                print(f"  {i}. Pose {r['pose_idx']}: RMSD={r['rmsd']:.3f} Å, GNINA={r['gnina_score']:.3f}")
            print("="*60)


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Compare docked poses to reference crystal structure")
    parser.add_argument('--receptor', required=True, help='Receptor PDB with reference ligand')
    parser.add_argument('--metadata', required=True, help='Receptor metadata JSON')
    parser.add_argument('--poses', required=True, help='Docked poses SDF file')
    parser.add_argument('--output', required=True, help='Output report path (without extension)')

    args = parser.parse_args()

    try:
        comparator = PoseComparator(
            receptor_pdb=args.receptor,
            metadata_path=args.metadata,
            docked_sdf=args.poses
        )

        results = comparator.compare_all_poses()
        comparator.generate_report(results, Path(args.output))

        print("Pose comparison complete!")

    except Exception as e:
        print(f"ERROR: Pose comparison failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
