#!/usr/bin/env python3
"""
Input Conformer Comparison
Compares input ligand conformer(s) to crystal structure ligand using alignment.
Tests if the starting conformer was close to the crystal pose.
"""

import sys
import json
from pathlib import Path
from typing import List, Dict, Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from Bio.PDB import PDBParser


class ConformerComparator:
    """Compares input conformers to reference crystal structure with alignment."""

    def __init__(self, receptor_pdb: Path, metadata_path: Path, input_sdf: Path):
        """
        Initialize comparator.

        Args:
            receptor_pdb: Path to aligned receptor PDB with reference ligand
            metadata_path: Path to receptor metadata JSON
            input_sdf: Path to input ligand SDF file
        """
        self.receptor_pdb = Path(receptor_pdb)
        self.metadata_path = Path(metadata_path)
        self.input_sdf = Path(input_sdf)

        # Load metadata
        with open(self.metadata_path) as f:
            self.metadata = json.load(f)

        self.ref_ligand_resn = self.metadata['target_ligand_resn']
        self.ref_ligand_chain = self.metadata['target_ligand_chain']

        print(f"Reference ligand: {self.ref_ligand_resn} chain {self.ref_ligand_chain}")

        # Extract reference ligand
        self.ref_mol = self._extract_reference_ligand()

        # Load input conformers
        self.input_mols = self._load_input_conformers()

        print(f"Loaded {len(self.input_mols)} input conformer(s)")

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

        # Store coords for later
        self.ref_coords = np.array(coords)
        self.ref_atoms = atoms

        # Load first input mol to use as template
        supplier = Chem.SDMolSupplier(str(self.input_sdf), removeHs=True, sanitize=True)
        template_mol = next(iter(supplier))

        if template_mol is None:
            raise ValueError("Could not load input ligand as template")

        if template_mol.GetNumAtoms() != len(coords):
            raise ValueError(f"Atom count mismatch: input={template_mol.GetNumAtoms()}, ref={len(coords)}")

        # Create reference mol with crystal coords
        ref_mol = Chem.Mol(template_mol)
        conf = ref_mol.GetConformer()
        for i in range(len(coords)):
            x, y, z = float(coords[i][0]), float(coords[i][1]), float(coords[i][2])
            conf.SetAtomPosition(i, (x, y, z))

        return ref_mol

    def _load_input_conformers(self) -> List[Chem.Mol]:
        """Load input conformer(s) from SDF file."""
        print(f"Loading input conformers from {self.input_sdf}")

        mols = []
        supplier = Chem.SDMolSupplier(str(self.input_sdf), removeHs=True, sanitize=True)

        for idx, mol in enumerate(supplier):
            if mol is None:
                print(f"WARNING: Could not parse conformer {idx}")
                continue
            mols.append(mol)

        return mols

    def compare_conformers(self) -> List[Dict]:
        """
        Calculate aligned RMSD for all input conformers vs reference.
        Uses RDKit's symmetry-aware alignment.

        Returns:
            List of result dicts sorted by RMSD
        """
        print("Calculating aligned RMSD for input conformers...")

        results = []

        for idx, mol in enumerate(self.input_mols):
            try:
                # Use symmetry-corrected alignment
                rmsd = rdMolAlign.GetBestRMS(mol, self.ref_mol)

                result = {
                    'conformer_idx': idx,
                    'rmsd_aligned': rmsd,
                    'pass_threshold': rmsd < 2.0
                }

                results.append(result)

                print(f"Conformer {idx}: Aligned RMSD={rmsd:.3f} Å")

            except Exception as e:
                print(f"ERROR: Failed to calculate RMSD for conformer {idx}: {e}")
                continue

        # Sort by RMSD
        results.sort(key=lambda x: x['rmsd_aligned'])

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
                'input_sdf': str(self.input_sdf),
                'total_conformers': len(results),
                'passing_conformers': sum(1 for r in results if r['pass_threshold']),
                'best_rmsd': results[0]['rmsd_aligned'] if results else None,
                'results': results
            }, f, indent=2)

        print(f"Saved JSON report to {json_path}")

        # CSV report
        csv_path = output_path.with_suffix('.csv')
        with open(csv_path, 'w') as f:
            f.write("conformer_idx,rmsd_aligned,pass_threshold\n")
            for r in results:
                f.write(f"{r['conformer_idx']},{r['rmsd_aligned']:.4f},{r['pass_threshold']}\n")

        print(f"Saved CSV report to {csv_path}")

        # Summary
        if results:
            print("\n" + "="*60)
            print("INPUT CONFORMER COMPARISON SUMMARY")
            print("="*60)
            print(f"Total conformers: {len(results)}")
            print(f"Passing conformers (RMSD < 2.0 Å): {sum(1 for r in results if r['pass_threshold'])}")
            print(f"\nBest conformer:")
            best = results[0]
            print(f"  Conformer {best['conformer_idx']}: Aligned RMSD={best['rmsd_aligned']:.3f} Å")
            if best['rmsd_aligned'] < 2.0:
                print(f"  ✓ Starting conformer was close to crystal pose")
            else:
                print(f"  ✗ Starting conformer was NOT close to crystal pose")
            print("="*60)


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Compare input ligand conformer to reference crystal structure")
    parser.add_argument('--receptor', required=True, help='Receptor PDB with reference ligand')
    parser.add_argument('--metadata', required=True, help='Receptor metadata JSON')
    parser.add_argument('--input', required=True, help='Input ligand SDF file')
    parser.add_argument('--output', required=True, help='Output report path (without extension)')

    args = parser.parse_args()

    try:
        comparator = ConformerComparator(
            receptor_pdb=args.receptor,
            metadata_path=args.metadata,
            input_sdf=args.input
        )

        results = comparator.compare_conformers()
        comparator.generate_report(results, Path(args.output))

        print("Conformer comparison complete!")

    except Exception as e:
        print(f"ERROR: Conformer comparison failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
