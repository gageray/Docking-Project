#!/usr/bin/env python3
"""
Input Conformer Comparison
Compares input ligand conformer(s) to crystal structure ligand using alignment.
Tests if the starting conformer was close to the crystal pose.
"""

import sys
import json
import zipfile
from pathlib import Path
from typing import List, Dict, Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


class ConformerComparator:
    """Compares input conformers to reference crystal structure with alignment."""

    def __init__(self, receptor_pdb: Path, metadata_path: Path, input_file: Path, smiles: str):
        """
        Initialize comparator.

        Args:
            receptor_pdb: Path to aligned receptor PDB with reference ligand
            metadata_path: Path to receptor metadata JSON
            input_file: Path to input ligand PDBQT zip file
            smiles: SMILES string for the ligand (to build template)
        """
        self.receptor_pdb = Path(receptor_pdb)
        self.metadata_path = Path(metadata_path)
        self.input_file = Path(input_file)
        self.smiles = smiles

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
        """Extract reference ligand from receptor PDB and map to SMILES topology."""
        print(f"Extracting reference ligand from {self.receptor_pdb}")

        # 1. Extract just the PDB lines for the reference ligand
        pdb_lines = []
        with open(self.receptor_pdb) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    resn = line[17:20].strip()
                    chain = line[21]
                    if resn == self.ref_ligand_resn and chain == self.ref_ligand_chain:
                        pdb_lines.append(line)

        if not pdb_lines:
            raise ValueError(f"Could not find ligand {self.ref_ligand_resn} in chain {self.ref_ligand_chain}")

        pdb_block = "".join(pdb_lines)

        # 2. Load the PDB block as a raw molecule (contains true 3D coords, but random atom order)
        ref_raw = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
        if ref_raw is None:
            raise ValueError("RDKit failed to parse the extracted PDB block.")

        print(f"Extracted {ref_raw.GetNumAtoms()} heavy atoms from reference ligand")

        # 3. Build the perfect template graph from SMILES
        template_mol = Chem.MolFromSmiles(self.smiles)
        if template_mol is None:
            raise ValueError(f"Invalid SMILES: {self.smiles}")

        # 4. The Magic Step: RDKit maps the SMILES graph onto the PDB coordinates!
        # This fixes the scrambled indices and guarantees isomorphic alignment.
        try:
            ref_mol = AllChem.AssignBondOrdersFromTemplate(template_mol, ref_raw)
        except Exception as e:
            raise ValueError(f"Could not map SMILES graph to PDB coordinates: {e}")

        return ref_mol

    def _load_input_conformers(self) -> List[Chem.Mol]:
        """Load input conformer(s) from PDBQT zip file."""
        print(f"Loading input conformers from {self.input_file}")

        # Build template from SMILES to get proper connectivity
        template = Chem.MolFromSmiles(self.smiles)
        if template is None:
            raise ValueError(f"Invalid SMILES: {self.smiles}")

        num_heavy = template.GetNumAtoms()
        mols = []

        with zipfile.ZipFile(self.input_file, 'r') as zf:
            pdbqt_files = [f for f in zf.namelist() if f.endswith('.pdbqt')]
            print(f"Found {len(pdbqt_files)} PDBQT files in zip")

            for pdbqt_name in sorted(pdbqt_files):
                pdbqt_content = zf.read(pdbqt_name).decode('utf-8')

                # 1. Strip AutoDock columns (67+) so RDKit can read it as standard PDB
                clean_pdb = []
                for line in pdbqt_content.splitlines():
                    if line.startswith(('ATOM', 'HETATM')):
                        clean_pdb.append(line[:66])
                    else:
                        clean_pdb.append(line)

                pdb_string = "\n".join(clean_pdb)

                # 2. Load the raw coordinates (with Meeko's scrambled atom order)
                raw_mol = Chem.MolFromPDBBlock(pdb_string, removeHs=True, sanitize=False)
                if raw_mol is None:
                    print(f"WARNING: RDKit failed to parse {pdbqt_name}")
                    continue

                if raw_mol.GetNumAtoms() != num_heavy:
                    print(f"WARNING: Atom count mismatch in {pdbqt_name}: {raw_mol.GetNumAtoms()} vs {num_heavy}")
                    continue

                # 3. Snap the perfect SMILES graph onto the Meeko coordinates
                try:
                    aligned_mol = AllChem.AssignBondOrdersFromTemplate(template, raw_mol)
                    mols.append(aligned_mol)
                except Exception as e:
                    print(f"WARNING: Failed to map SMILES to {pdbqt_name}: {e}")

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
                'input_file': str(self.input_file),
                'smiles': self.smiles,
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
    parser.add_argument('--input', required=True, help='Input ligand PDBQT zip file')
    parser.add_argument('--smiles', required=True, help='SMILES string for ligand')
    parser.add_argument('--output', required=True, help='Output report path (without extension)')

    args = parser.parse_args()

    try:
        comparator = ConformerComparator(
            receptor_pdb=args.receptor,
            metadata_path=args.metadata,
            input_file=args.input,
            smiles=args.smiles
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
