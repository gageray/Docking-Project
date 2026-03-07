#!/usr/bin/env python3
"""
CREST-based Conformer Generation Validation
Uses GFN2-xTB via CREST for accurate conformer generation of charged species.

Workflow:
1. Generate neutral seed with ETKDG (bypass charged embedding failure)
2. Map neutral heavy atoms to charged graph to safely add proton in 3D space
3. Export to XYZ for CREST
4. Dynamic CREST execution: route to `xtb --opt` for rigid, or `crest` for flexible
5. Parse ensemble and compare to crystal structure
"""

import sys
import json
import subprocess
import tempfile
import shutil
import os
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, Descriptors, rdMolDescriptors, rdFMCS
from Bio.PDB import PDBParser


class CRESTConformerValidator:
    """Validates conformer generation using CREST/GFN2-xTB."""

    def __init__(self, neutral_smiles: str, charged_smiles: str, energy_window: float = 5.0, nocross: bool = False):
        """
        Initialize validator.

        Args:
            neutral_smiles: Neutral SMILES for seed generation
            charged_smiles: Charged SMILES (pH 7.4) for final structures
            energy_window: Energy window in kcal/mol (CREST uses --ewin)
            nocross: Disable genetic crossing for rigid molecules
        """
        self.neutral_smiles = neutral_smiles
        self.charged_smiles = charged_smiles
        self.energy_window = energy_window
        self.nocross = nocross

        # Check CREST availability
        if not shutil.which('crest') and not shutil.which('xtb'):
            raise RuntimeError("CREST or xTB not found in PATH.")

        # Get formal charge from charged SMILES
        mol_charged = Chem.MolFromSmiles(charged_smiles)
        self.formal_charge = Chem.GetFormalCharge(mol_charged)

        print(f"Neutral SMILES: {neutral_smiles}")
        print(f"Charged SMILES: {charged_smiles}")
        print(f"Formal charge: {self.formal_charge:+d}")

    def generate_seed(self, output_xyz: Path) -> Chem.Mol:
        """Generate 3D charged seed by embedding neutral graph and transferring coords."""
        print("\nGenerating 3D seed via neutral embedding...")

        # 1. Embed Neutral to bypass ETKDG charge failures
        neutral_mol = Chem.MolFromSmiles(self.neutral_smiles)
        if neutral_mol is None:
            raise ValueError(f"Invalid neutral SMILES: {self.neutral_smiles}")
        
        neutral_mol = Chem.AddHs(neutral_mol)
        
        # Single conformer with ETKDG on the safe neutral geometry
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(neutral_mol, params) != 0:
            raise RuntimeError("ETKDG embedding failed for neutral molecule")
        
        # MMFF optimization on the neutral molecule
        if AllChem.MMFFOptimizeMolecule(neutral_mol) != 0:
            print("WARNING: MMFF optimization did not converge")

        # 2. Map neutral heavy atoms to charged heavy atoms
        charged_mol = Chem.MolFromSmiles(self.charged_smiles)
        if charged_mol is None:
            raise ValueError(f"Invalid charged SMILES: {self.charged_smiles}")

        neutral_heavy = Chem.RemoveHs(neutral_mol)
        charged_heavy = Chem.RemoveHs(charged_mol)

        # Ignore charges and bond order shifts to find the matching core using MCS
        mcs = rdFMCS.FindMCS(
            [neutral_heavy, charged_heavy],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            matchValences=False
        )
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        match_neutral = neutral_heavy.GetSubstructMatch(mcs_mol)
        match_charged = charged_heavy.GetSubstructMatch(mcs_mol)

        if len(match_neutral) != neutral_heavy.GetNumAtoms():
            raise RuntimeError("Could not map neutral and charged heavy atoms.")

        # 3. Transfer 3D coordinates from the neutral graph to the charged graph
        conf = Chem.Conformer(charged_heavy.GetNumAtoms())
        neutral_conf = neutral_heavy.GetConformer()
        for n_idx, c_idx in zip(match_neutral, match_charged):
            conf.SetAtomPosition(c_idx, neutral_conf.GetAtomPosition(n_idx))
        charged_heavy.AddConformer(conf)

        # 4. Add Hydrogens (including the new proton) mathematically in 3D space
        charged_3d = Chem.AddHs(charged_heavy, addCoords=True)
        
        # Relax the newly added proton into a local minimum
        AllChem.MMFFOptimizeMolecule(charged_3d)

        # 5. Export the correct fully-protonated (e.g., 37-atom) charged seed to XYZ for CREST
        conf_3d = charged_3d.GetConformer()
        with open(output_xyz, 'w') as f:
            f.write(f"{charged_3d.GetNumAtoms()}\n")
            f.write(f"Charged seed for CREST\n")
            for i in range(charged_3d.GetNumAtoms()):
                atom = charged_3d.GetAtomWithIdx(i)
                pos = conf_3d.GetAtomPosition(i)
                f.write(f"{atom.GetSymbol():2s} {pos.x:12.6f} {pos.y:12.6f} {pos.z:12.6f}\n")

        print(f"Wrote {charged_3d.GetNumAtoms()}-atom charged seed to {output_xyz}")
        return charged_3d

    def run_crest(self, seed_xyz: Path, output_dir: Path, num_threads: int, n_rot: int, search_level: str = 'quick') -> Tuple[Path, List[float]]:
        """Run CREST conformational search or xTB optimization dynamically based on flexibility."""
        print(f"\nRunning GFN2-xTB calculation...")
        print(f"  Formal charge: {self.formal_charge:+d}")
        print(f"  Rotatable bonds: {n_rot}")
        print(f"  Energy window: {self.energy_window} kcal/mol")
        print(f"  Threads: {num_threads}")

        # Dynamic routing to prevent metadynamics bias trap on rigid molecules
        if n_rot == 0:
            # Completely rigid: just relax the charged state, no conformational search needed
            cmd = [
                'xtb',
                str(seed_xyz),
                '--opt',
                '--chrg', str(self.formal_charge),
                '-T', str(num_threads)
            ]
        elif n_rot <= 2:
            # Stiff but rotatable: Fast search, append nocross if flagged
            cmd = [
                'crest',
                str(seed_xyz),
                '--gfn2',
                '--chrg', str(self.formal_charge),
                '--ewin', str(self.energy_window)
            ]
            if self.nocross:
                cmd.append('--nocross')
            if search_level in ['quick', 'squick', 'mquick']:
                cmd.append(f'--{search_level}')
            cmd.extend(['-T', str(num_threads)])
        else:
            # Flexible molecule: Standard quick search
            cmd = [
                'crest',
                str(seed_xyz),
                '--gfn2',
                '--chrg', str(self.formal_charge),
                '--ewin', str(self.energy_window)
            ]
            if search_level in ['quick', 'squick', 'mquick']:
                cmd.append(f'--{search_level}')
            cmd.extend(['-T', str(num_threads)])

        print(f"  Command: {' '.join(cmd)}")

        try:
            # Setup environment to prevent underlying math libraries from conflicting with CREST's threading
            crest_env = os.environ.copy()
            crest_env['OMP_NUM_THREADS'] = '1'
            crest_env['MKL_NUM_THREADS'] = '1'
            crest_env['OPENBLAS_NUM_THREADS'] = '1'

            # Run calculation and stream output to console
            process = subprocess.Popen(
                cmd,
                cwd=output_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                env=crest_env
            )
            
            # Print output as it comes in to monitor SCF iterations
            for line in process.stdout:
                print(line, end='')
                
            process.wait()

            if process.returncode != 0:
                raise RuntimeError(f"Execution failed with return code {process.returncode}")

            print("\nCalculation completed successfully")

        except Exception as e:
            raise RuntimeError(f"Execution failed: {e}")

        # Parse output files
        conformers_file = output_dir / "crest_conformers.xyz"
        energies_file = output_dir / "crest.energies"

        # If we ran xtb --opt instead of crest, mock the crest outputs so the parser doesn't crash
        if n_rot == 0:
            xtb_out = output_dir / "xtbopt.xyz"
            if xtb_out.exists():
                shutil.copy(xtb_out, conformers_file)
                with open(energies_file, 'w') as f:
                    f.write("0.00000\n")

        if not conformers_file.exists():
            raise RuntimeError(f"Output not found: {conformers_file}")

        # Read energies
        energies = []
        if energies_file.exists():
            with open(energies_file) as f:
                for line in f:
                    # Energy file format: relative_energy(kcal/mol)
                    parts = line.strip().split()
                    if parts:
                        energies.append(float(parts[0]))

        print(f"Found {len(energies)} conformers")
        return conformers_file, energies

    def parse_crest_ensemble(self, xyz_file: Path, energies: List[float]) -> List[Chem.Mol]:
        """Parse CREST multi-XYZ file and assign topology from charged SMILES."""
        print(f"\nParsing CREST ensemble from {xyz_file}...")

        # Template molecule from charged SMILES
        template = Chem.MolFromSmiles(self.charged_smiles)
        if template is None:
            raise ValueError(f"Invalid charged SMILES: {self.charged_smiles}")

        template = Chem.AddHs(template)
        num_atoms = template.GetNumAtoms()

        conformers = []
        with open(xyz_file) as f:
            while True:
                # Read header
                natoms_line = f.readline()
                if not natoms_line:
                    break

                natoms = int(natoms_line.strip())
                comment = f.readline()

                if natoms != num_atoms:
                    raise ValueError(f"Atom count mismatch: XYZ has {natoms}, template has {num_atoms}")

                # Read coordinates
                coords = []
                for _ in range(natoms):
                    line = f.readline()
                    parts = line.strip().split()
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                    coords.append([x, y, z])

                # Create mol with template topology
                mol = Chem.Mol(template)
                conf = Chem.Conformer(num_atoms)
                for i, (x, y, z) in enumerate(coords):
                    conf.SetAtomPosition(i, (x, y, z))
                mol.AddConformer(conf, assignId=True)

                conformers.append(mol)

        print(f"Parsed {len(conformers)} conformers")
        return conformers

    def compare_to_crystal(self, conformers: List[Chem.Mol], energies: List[float],
                          receptor_pdb: Path, metadata_path: Path) -> Dict:
        """Compare conformers to crystal structure."""
        print("\nComparing to crystal structure...")

        # Load metadata
        with open(metadata_path) as f:
            metadata = json.load(f)

        ref_ligand_resn = metadata['target_ligand_resn']
        ref_ligand_chain = metadata['target_ligand_chain']

        # Extract reference ligand from PDB
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

        # Extract coords (heavy atoms only)
        coords = []
        for atom in ligand_residue:
            if atom.element.strip() != 'H':
                coords.append(atom.coord)

        # Create reference mol using the generated conformer as a topological template
        ref_mol = Chem.RemoveHs(conformers[0])
        ref_conf = ref_mol.GetConformer(0)

        if len(coords) != ref_mol.GetNumAtoms():
            raise ValueError(f"Atom count mismatch: generated={ref_mol.GetNumAtoms()}, crystal={len(coords)}")

        for i in range(len(coords)):
            x, y, z = float(coords[i][0]), float(coords[i][1]), float(coords[i][2])
            ref_conf.SetAtomPosition(i, (x, y, z))

        # Compare all conformers (strip H for comparison)
        results = []
        min_energy = min(energies) if energies else 0.0

        for idx, (mol, energy) in enumerate(zip(conformers, energies)):
            mol_no_h = Chem.RemoveHs(mol)
            rmsd = rdMolAlign.GetBestRMS(mol_no_h, ref_mol, prbId=0, refId=0)

            delta_e = energy - min_energy

            results.append({
                'conf_id': idx,
                'energy': float(energy),
                'delta_energy': float(delta_e),
                'rmsd_to_crystal': float(rmsd),
                'pass_threshold': bool(rmsd < 2.0)
            })

        # Sort by RMSD
        results.sort(key=lambda x: x['rmsd_to_crystal'])

        best = results[0]
        summary = {
            'reference_ligand': f"{ref_ligand_resn}:{ref_ligand_chain}",
            'method': 'CREST/GFN2-xTB',
            'total_conformers': len(conformers),
            'best_rmsd': best['rmsd_to_crystal'],
            'best_conf_id': best['conf_id'],
            'best_conf_delta_energy': best['delta_energy'],
            'passing_conformers': sum(1 for r in results if r['pass_threshold']),
            'all_results': results
        }

        print(f"Best conformer: ID={best['conf_id']}, RMSD={best['rmsd_to_crystal']:.3f} Å")
        print(f"Energy: ΔE={best['delta_energy']:.2f} kcal/mol")
        print(f"Conformers within 2.0 Å: {summary['passing_conformers']}/{len(results)}")

        if best['rmsd_to_crystal'] < 2.0:
            print("✓ Successfully recovered crystal-like pose")
        else:
            print("✗ Failed to recover crystal-like pose")

        return summary

    def save_conformers(self, conformers: List[Chem.Mol], energies: List[float], output_sdf: Path):
        """Save conformers to SDF."""
        writer = Chem.SDWriter(str(output_sdf))

        min_energy = min(energies) if energies else 0.0

        for idx, (mol, energy) in enumerate(zip(conformers, energies)):
            mol.SetProp('Energy', str(energy))
            mol.SetProp('DeltaEnergy', str(energy - min_energy))
            mol.SetProp('ConformerID', str(idx))
            mol.SetProp('Method', 'CREST/GFN2-xTB')
            writer.write(mol, confId=0)

        writer.close()
        print(f"\nSaved {len(conformers)} conformers to {output_sdf}")


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Validate conformer generation using CREST/GFN2-xTB")
    parser.add_argument('--neutral-smiles', required=True, help='Neutral SMILES for seed generation')
    parser.add_argument('--charged-smiles', required=True, help='Charged SMILES (pH 7.4)')
    parser.add_argument('--energy-window', type=float, default=5.0, help='Energy window in kcal/mol')
    parser.add_argument('--search-level', type=str, choices=['default', 'quick', 'squick', 'mquick'], default='quick', help='CREST search thoroughness')
    parser.add_argument('--nocross', action='store_true', help='Disable genetic crossing (recommended for rigid molecules)')
    parser.add_argument('--receptor', required=True, help='Receptor PDB with reference ligand')
    parser.add_argument('--metadata', required=True, help='Receptor metadata JSON')
    parser.add_argument('--output', required=True, help='Output report path (without extension)')
    parser.add_argument('--save-sdf', required=True, help='Save conformers to SDF file')
    parser.add_argument('--threads', type=int, default=8, help='Number of CPU threads for CREST')
    parser.add_argument('--keep-crest-files', action='store_true', help='Keep CREST working directory')

    args = parser.parse_args()

    try:
        validator = CRESTConformerValidator(
            neutral_smiles=args.neutral_smiles,
            charged_smiles=args.charged_smiles,
            energy_window=args.energy_window,
            nocross=args.nocross
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            seed_xyz = tmpdir / "seed.xyz"

            # Generate seed via neutral embedding and proton mapping
            mol_3d = validator.generate_seed(seed_xyz)
            
            # Calculate rotatable bonds to dynamically route crest/xtb execution
            n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol_3d)

            # Run CREST or xTB
            conformers_file, energies = validator.run_crest(seed_xyz, tmpdir, args.threads, n_rot, args.search_level)

            # Parse ensemble
            conformers = validator.parse_crest_ensemble(conformers_file, energies)

            # Compare to crystal
            results = validator.compare_to_crystal(
                conformers,
                energies,
                Path(args.receptor),
                Path(args.metadata)
            )

            # Save conformers
            validator.save_conformers(conformers, energies, Path(args.save_sdf))

            # Save report
            output_path = Path(args.output)
            json_path = output_path.with_suffix('.json')
            with open(json_path, 'w') as f:
                json.dump({
                    'neutral_smiles': args.neutral_smiles,
                    'charged_smiles': args.charged_smiles,
                    'formal_charge': validator.formal_charge,
                    'energy_window': args.energy_window,
                    'method': 'CREST/GFN2-xTB',
                    **results
                }, f, indent=2)

            print(f"\nSaved report to {json_path}")

            if args.keep_crest_files:
                crest_dir = output_path.parent / f"{output_path.stem}_crest_files"
                shutil.copytree(tmpdir, crest_dir)
                print(f"Saved CREST files to {crest_dir}")

        print("\nValidation complete!")

    except Exception as e:
        print(f"ERROR: Validation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()