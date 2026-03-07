#!/usr/bin/env python3
"""
Conformer Generation Validation
Validates that conformer generation from SMILES produces reasonable, diverse conformers
and tests if any match the crystal structure pose (for known ligands).

Checks:
1. Geometric validity (no steric clashes, valid bond geometry)
2. Energetic reasonableness (MMFF energy within thermal window)
3. Structural diversity (RMSD matrix between conformers)
4. Crystal pose recovery (if reference structure provided)
"""

import sys
import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, Descriptors, rdForceFieldHelpers
from Bio.PDB import PDBParser


class ConformerValidator:
    """Validates conformer generation quality."""

    def __init__(self, smiles: str, num_conformers: int = 50, energy_window: float = 40.0):
        """
        Initialize validator.

        Args:
            smiles: SMILES string
            num_conformers: Number of conformers to generate
            energy_window: Energy window in kcal/mol above minimum (default 40.0 for induced-fit)
        """
        self.smiles = smiles
        self.num_conformers = num_conformers
        self.energy_window = energy_window

        print(f"Generating {num_conformers} conformers from SMILES: {smiles}")

        # Generate conformers
        self.mol = self._generate_conformers()

        # Calculate energies
        self.energies = self._calculate_energies()

        # Store min energy for reference
        self.min_energy = min(e for e in self.energies if e != float('inf'))

        # Filter by energy
        self.valid_conf_ids = self._filter_by_energy()

        # Print energy distribution
        self._print_energy_distribution()

        print(f"Generated {self.mol.GetNumConformers()} conformers")
        print(f"Min energy: {self.min_energy:.2f} kcal/mol")
        print(f"Valid conformers (ΔE < {energy_window} kcal/mol): {len(self.valid_conf_ids)}")

    def _generate_conformers(self) -> Chem.Mol:
        """Generate conformers using ETKDG."""
        mol = Chem.MolFromSmiles(self.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {self.smiles}")

        mol = Chem.AddHs(mol)

        # ETKDG conformer generation
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.numThreads = 0

        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=self.num_conformers,
            params=params
        )

        if len(conf_ids) == 0:
            raise ValueError("Failed to generate any conformers")

        # MMFF optimization
        print("Optimizing conformers with MMFF...")
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)

        # Check for failed optimizations
        failed = sum(1 for converged, energy in results if converged != 0)
        if failed > 0:
            print(f"WARNING: {failed} conformers failed MMFF optimization")

        # Remove hydrogens for comparison
        mol = Chem.RemoveHs(mol)

        return mol

    def _calculate_energies(self) -> List[float]:
        """Calculate MMFF energies for all conformers."""
        energies = []

        # Need H for energy calculation
        mol_h = Chem.AddHs(self.mol)

        for conf_id in range(self.mol.GetNumConformers()):
            # Copy conformer to mol with H
            conf = self.mol.GetConformer(conf_id)
            conf_h = mol_h.GetConformer(conf_id)

            # Get MMFF energy
            props = AllChem.MMFFGetMoleculeProperties(mol_h)
            ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=conf_id)

            if ff is None:
                print(f"WARNING: Could not get force field for conformer {conf_id}")
                energies.append(float('inf'))
            else:
                energy = ff.CalcEnergy()
                energies.append(energy)

        return energies

    def _filter_by_energy(self) -> List[int]:
        """Filter conformers by energy window."""
        min_energy = min(e for e in self.energies if e != float('inf'))

        valid_ids = []
        for conf_id, energy in enumerate(self.energies):
            if energy != float('inf') and (energy - min_energy) <= self.energy_window:
                valid_ids.append(conf_id)

        return valid_ids

    def _print_energy_distribution(self):
        """Print energy distribution in bins."""
        print("\nEnergy Distribution:")

        # Create 10 kcal/mol bins
        bins = [0, 10, 20, 30, 40, 50, 100, float('inf')]
        bin_counts = {f"{bins[i]}-{bins[i+1]}": 0 for i in range(len(bins)-1)}

        for energy in self.energies:
            if energy == float('inf'):
                bin_counts[f"{bins[-2]}-inf"] += 1
                continue

            delta_e = energy - self.min_energy
            for i in range(len(bins)-1):
                if delta_e >= bins[i] and delta_e < bins[i+1]:
                    bin_counts[f"{bins[i]}-{bins[i+1]}"] += 1
                    break

        for bin_range, count in bin_counts.items():
            if count > 0:
                print(f"  ΔE {bin_range} kcal/mol: {count} conformers")

    def check_geometry(self) -> Dict:
        """Check geometric validity of conformers."""
        print("\nChecking geometric validity...")

        results = {
            'num_atoms': self.mol.GetNumAtoms(),
            'num_bonds': self.mol.GetNumBonds(),
            'num_rotatable_bonds': Descriptors.NumRotatableBonds(self.mol),
            'molecular_weight': Descriptors.MolWt(self.mol),
            'conformers': []
        }

        for conf_id in self.valid_conf_ids:
            conf = self.mol.GetConformer(conf_id)

            # Check for steric clashes (atoms too close)
            coords = conf.GetPositions()
            min_dist = float('inf')

            for i in range(len(coords)):
                for j in range(i+1, len(coords)):
                    # Skip bonded atoms
                    if self.mol.GetBondBetweenAtoms(i, j) is not None:
                        continue

                    dist = np.linalg.norm(coords[i] - coords[j])
                    min_dist = min(min_dist, dist)

            # Flag clashes (atoms closer than 1.5 Å)
            has_clash = min_dist < 1.5

            results['conformers'].append({
                'conf_id': int(conf_id),
                'energy': float(self.energies[conf_id]),
                'min_nonbonded_dist': float(min_dist),
                'has_clash': bool(has_clash)
            })

        num_clashes = sum(1 for c in results['conformers'] if c['has_clash'])
        print(f"Conformers with steric clashes: {num_clashes}/{len(self.valid_conf_ids)}")

        return results

    def check_diversity(self) -> Dict:
        """Calculate RMSD matrix to check conformational diversity."""
        print("\nChecking conformational diversity...")

        n_conf = len(self.valid_conf_ids)

        # Handle single conformer case
        if n_conf < 2:
            print("WARNING: Only 1 conformer, cannot assess diversity")
            return {
                'num_conformers': n_conf,
                'mean_pairwise_rmsd': None,
                'median_pairwise_rmsd': None,
                'min_pairwise_rmsd': None,
                'max_pairwise_rmsd': None,
                'rmsd_matrix': []
            }

        rmsd_matrix = np.zeros((n_conf, n_conf))

        for i, conf_id_i in enumerate(self.valid_conf_ids):
            for j, conf_id_j in enumerate(self.valid_conf_ids):
                if i < j:
                    rmsd = rdMolAlign.GetBestRMS(
                        self.mol, self.mol,
                        prbId=conf_id_i,
                        refId=conf_id_j
                    )
                    rmsd_matrix[i, j] = rmsd
                    rmsd_matrix[j, i] = rmsd

        # Statistics
        upper_triangle = rmsd_matrix[np.triu_indices_from(rmsd_matrix, k=1)]

        results = {
            'num_conformers': n_conf,
            'mean_pairwise_rmsd': float(np.mean(upper_triangle)),
            'median_pairwise_rmsd': float(np.median(upper_triangle)),
            'min_pairwise_rmsd': float(np.min(upper_triangle)),
            'max_pairwise_rmsd': float(np.max(upper_triangle)),
            'rmsd_matrix': rmsd_matrix.tolist()
        }

        if results['mean_pairwise_rmsd'] is not None:
            print(f"Mean pairwise RMSD: {results['mean_pairwise_rmsd']:.3f} Å")
            print(f"RMSD range: {results['min_pairwise_rmsd']:.3f} - {results['max_pairwise_rmsd']:.3f} Å")

            # Check if diverse (mean RMSD > 1.0 Å suggests good diversity)
            if results['mean_pairwise_rmsd'] < 1.0:
                print("WARNING: Low conformational diversity (mean RMSD < 1.0 Å)")
            else:
                print("✓ Good conformational diversity")

        return results

    def compare_to_crystal(self, receptor_pdb: Path, metadata_path: Path) -> Dict:
        """Compare conformers to crystal structure ligand."""
        print("\nComparing to crystal structure...")

        # Load metadata
        with open(metadata_path) as f:
            metadata = json.load(f)

        ref_ligand_resn = metadata['target_ligand_resn']
        ref_ligand_chain = metadata['target_ligand_chain']

        # Extract reference ligand
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('receptor', str(receptor_pdb))

        ligand_residue = None
        for model in structure:
            for chain in model:
                if chain.id == ref_ligand_chain:
                    for residue in chain:
                        if residue.resname == ref_ligand_resn:
                            ligand_residue = residue
                            break

        if not ligand_residue:
            raise ValueError(f"Could not find ligand {ref_ligand_resn}")

        # Extract coords
        coords = []
        for atom in ligand_residue:
            if atom.element.strip() != 'H':
                coords.append(atom.coord)

        # Create reference mol with crystal coords
        ref_mol = Chem.Mol(self.mol)
        ref_conf = ref_mol.GetConformer(0)

        if len(coords) != self.mol.GetNumAtoms():
            raise ValueError(f"Atom count mismatch: generated={self.mol.GetNumAtoms()}, crystal={len(coords)}")

        for i in range(len(coords)):
            x, y, z = float(coords[i][0]), float(coords[i][1]), float(coords[i][2])
            ref_conf.SetAtomPosition(i, (x, y, z))

        # Compare all conformers to crystal
        results = []
        for conf_id in self.valid_conf_ids:
            rmsd = rdMolAlign.GetBestRMS(self.mol, ref_mol, prbId=conf_id, refId=0)

            energy = self.energies[conf_id]
            delta_e = energy - self.min_energy

            # Determine energy bin
            if delta_e < 10:
                energy_bin = "0-10"
            elif delta_e < 20:
                energy_bin = "10-20"
            elif delta_e < 30:
                energy_bin = "20-30"
            elif delta_e < 40:
                energy_bin = "30-40"
            else:
                energy_bin = "40+"

            results.append({
                'conf_id': int(conf_id),
                'energy': float(energy),
                'delta_energy': float(delta_e),
                'energy_bin': energy_bin,
                'rmsd_to_crystal': float(rmsd),
                'pass_threshold': bool(rmsd < 2.0)
            })

        # Sort by RMSD
        results.sort(key=lambda x: x['rmsd_to_crystal'])

        # Find energy rank of best RMSD conformer
        best_conf_id = results[0]['conf_id']
        energy_sorted = sorted(enumerate(self.energies), key=lambda x: x[1])
        energy_rank = next(i for i, (conf_id, _) in enumerate(energy_sorted) if conf_id == best_conf_id)
        best_delta_e = results[0]['delta_energy']

        summary = {
            'reference_ligand': f"{ref_ligand_resn}:{ref_ligand_chain}",
            'best_rmsd': results[0]['rmsd_to_crystal'],
            'best_conf_id': best_conf_id,
            'best_conf_energy_rank': energy_rank,
            'best_conf_delta_energy': best_delta_e,
            'best_conf_energy_bin': results[0]['energy_bin'],
            'passing_conformers': sum(1 for r in results if r['pass_threshold']),
            'all_results': results
        }

        print(f"Best conformer: ID={best_conf_id}, RMSD={results[0]['rmsd_to_crystal']:.3f} Å")
        print(f"Energy: ΔE={best_delta_e:.2f} kcal/mol (bin: {results[0]['energy_bin']} kcal/mol)")
        print(f"Energy rank: {energy_rank}/{len(self.energies)}")
        print(f"Conformers within 2.0 Å of crystal: {summary['passing_conformers']}/{len(results)}")

        if results[0]['rmsd_to_crystal'] < 2.0:
            print("✓ Successfully recovered crystal-like pose")
        else:
            print("✗ Failed to recover crystal-like pose")

        # Analyze if induced-fit binding
        if best_delta_e > 20.0:
            print(f"\n⚠ INDUCED-FIT BINDING DETECTED:")
            print(f"  Crystal pose is {best_delta_e:.1f} kcal/mol above solution minimum")
            print(f"  Binding site must provide ≥{best_delta_e:.1f} kcal/mol stabilization")

        return summary

    def generate_report(self, geometry_results: Dict, diversity_results: Dict,
                       crystal_results: Optional[Dict], output_path: Path):
        """Generate validation report."""
        output_path = Path(output_path)

        report = {
            'smiles': self.smiles,
            'num_conformers_generated': self.mol.GetNumConformers(),
            'num_conformers_valid': len(self.valid_conf_ids),
            'energy_window': self.energy_window,
            'min_energy': min(self.energies),
            'geometry': geometry_results,
            'diversity': diversity_results,
        }

        if crystal_results:
            report['crystal_comparison'] = crystal_results

        # JSON report
        json_path = output_path.with_suffix('.json')
        with open(json_path, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"\nSaved JSON report to {json_path}")

        # CSV for conformer data
        csv_path = output_path.with_suffix('.csv')
        with open(csv_path, 'w') as f:
            headers = ['conf_id', 'energy', 'delta_energy', 'energy_bin', 'min_nonbonded_dist', 'has_clash']
            if crystal_results:
                headers.extend(['rmsd_to_crystal', 'pass_threshold'])
            f.write(','.join(headers) + '\n')

            min_energy = min(self.energies)
            for i, conf_id in enumerate(self.valid_conf_ids):
                energy = self.energies[conf_id]
                delta_e = energy - min_energy
                geom = geometry_results['conformers'][i]

                # Determine energy bin
                if delta_e < 10:
                    energy_bin = "0-10"
                elif delta_e < 20:
                    energy_bin = "10-20"
                elif delta_e < 30:
                    energy_bin = "20-30"
                elif delta_e < 40:
                    energy_bin = "30-40"
                else:
                    energy_bin = "40+"

                row = [
                    str(conf_id),
                    f"{energy:.4f}",
                    f"{delta_e:.4f}",
                    energy_bin,
                    f"{geom['min_nonbonded_dist']:.4f}",
                    str(geom['has_clash'])
                ]

                if crystal_results:
                    crystal_data = next(r for r in crystal_results['all_results'] if r['conf_id'] == conf_id)
                    row.append(f"{crystal_data['rmsd_to_crystal']:.4f}")
                    row.append(str(crystal_data['pass_threshold']))

                f.write(','.join(row) + '\n')

        print(f"Saved CSV report to {csv_path}")

        # Summary
        print("\n" + "="*60)
        print("CONFORMER GENERATION VALIDATION SUMMARY")
        print("="*60)
        print(f"SMILES: {self.smiles}")
        print(f"Conformers generated: {self.mol.GetNumConformers()}")
        print(f"Valid conformers: {len(self.valid_conf_ids)}")
        print(f"Rotatable bonds: {geometry_results['num_rotatable_bonds']}")
        print(f"\nGeometry:")
        print(f"  Steric clashes: {sum(1 for c in geometry_results['conformers'] if c['has_clash'])}/{len(self.valid_conf_ids)}")
        print(f"\nDiversity:")
        if diversity_results['mean_pairwise_rmsd'] is not None:
            print(f"  Mean pairwise RMSD: {diversity_results['mean_pairwise_rmsd']:.3f} Å")
        else:
            print(f"  Mean pairwise RMSD: N/A (only 1 conformer)")
        if crystal_results:
            print(f"\nCrystal comparison:")
            print(f"  Best RMSD: {crystal_results['best_rmsd']:.3f} Å")
            print(f"  Energy rank: {crystal_results['best_conf_energy_rank']}/{len(self.energies)}")
            print(f"  Passing conformers: {crystal_results['passing_conformers']}/{len(self.valid_conf_ids)}")
        print("="*60)

    def save_conformers(self, output_sdf: Path):
        """Save valid conformers to SDF file."""
        writer = Chem.SDWriter(str(output_sdf))

        min_energy = min(self.energies)
        for conf_id in self.valid_conf_ids:
            # Write conformer with energy data
            writer.write(self.mol, confId=conf_id)
            mol_copy = Chem.Mol(self.mol, confId=conf_id)
            mol_copy.SetProp('Energy', str(self.energies[conf_id]))
            mol_copy.SetProp('DeltaEnergy', str(self.energies[conf_id] - min_energy))
            mol_copy.SetProp('ConformerID', str(conf_id))

        writer.close()
        print(f"\nSaved {len(self.valid_conf_ids)} conformers to {output_sdf}")


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Validate conformer generation from SMILES")
    parser.add_argument('--smiles', required=True, help='SMILES string')
    parser.add_argument('--num-conformers', type=int, default=50, help='Number of conformers to generate')
    parser.add_argument('--energy-window', type=float, default=20.0, help='Energy window in kcal/mol')
    parser.add_argument('--receptor', help='Optional: Receptor PDB with reference ligand for comparison')
    parser.add_argument('--metadata', help='Optional: Receptor metadata JSON')
    parser.add_argument('--output', required=True, help='Output report path (without extension)')
    parser.add_argument('--save-sdf', help='Optional: Save conformers to SDF file')

    args = parser.parse_args()

    try:
        validator = ConformerValidator(
            smiles=args.smiles,
            num_conformers=args.num_conformers,
            energy_window=args.energy_window
        )

        geometry_results = validator.check_geometry()
        diversity_results = validator.check_diversity()

        crystal_results = None
        if args.receptor and args.metadata:
            crystal_results = validator.compare_to_crystal(
                receptor_pdb=Path(args.receptor),
                metadata_path=Path(args.metadata)
            )

        validator.generate_report(
            geometry_results,
            diversity_results,
            crystal_results,
            Path(args.output)
        )

        if args.save_sdf:
            validator.save_conformers(Path(args.save_sdf))

        print("\nValidation complete!")

    except Exception as e:
        print(f"ERROR: Validation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
