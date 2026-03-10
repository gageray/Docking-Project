#!/usr/bin/env python3
"""
Pose Comparison and Validation
Compares docked poses against reference crystal structure ligand using RMSD and interaction analysis.
"""

import sys
import os
import json
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from Bio.PDB import PDBParser, Selection

# Add kaggle scripts to path for Google Drive auth
sys.path.insert(0, str(Path(__file__).parent.parent / 'kaggle'))
try:
    import drive_auth
    import drive_io
except ImportError:
    print("Warning: Could not import Google Drive modules")
    drive_auth = None
    drive_io = None


class PoseComparator:
    """Compares docked poses to reference crystal structure."""

    def __init__(self, receptor_pdb: Path, metadata_path: Path, docked_path: Path):
        """
        Initialize comparator.

        Args:
            receptor_pdb: Path to aligned receptor PDB with reference ligand
            metadata_path: Path to receptor metadata JSON
            docked_path: Path to GNINA output SDF with docked poses, or directory of SDFs
        """
        self.receptor_pdb = Path(receptor_pdb)
        self.metadata_path = Path(metadata_path)
        self.docked_path = Path(docked_path)

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
        """Load all docked poses from SDF file or directory of SDFs."""
        print(f"Loading docked poses from {self.docked_path}")

        poses = []
        sdf_files = []
        
        if self.docked_path.is_file():
            sdf_files.append(self.docked_path)
        elif self.docked_path.is_dir():
            sdf_files.extend(list(self.docked_path.rglob("*.sdf")))
            
        if not sdf_files:
            print(f"WARNING: No SDF files found in {self.docked_path}")
            return poses

        for sdf_file in sdf_files:
            supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=True, sanitize=True)

            for idx, mol in enumerate(supplier):
                if mol is None:
                    print(f"WARNING: Could not parse molecule {idx} in {sdf_file.name}")
                    continue

                # Extract GNINA properties
                props = {
                    'source_file': sdf_file.name,
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

        IMPORTANT: For docking validation, this uses RAW (unaligned) RMSD
        because both molecules should already be in the same coordinate frame.
        Alignment is only appropriate for conformer comparison, not pose validation.

        Args:
            probe_mol: Docked pose molecule
            ref_mol: Reference molecule (if None, uses self.ref_coords)

        Returns:
            RMSD in Angstroms (raw, no alignment)
        """
        # Extract probe coordinates
        probe_conf = probe_mol.GetConformer()
        probe_coords = np.array([probe_conf.GetAtomPosition(i) for i in range(probe_mol.GetNumAtoms())])

        # Get reference coordinates
        if ref_mol is not None:
            ref_conf = ref_mol.GetConformer()
            ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in range(ref_mol.GetNumAtoms())])
        else:
            ref_coords = self.ref_coords

        if probe_coords.shape[0] != ref_coords.shape[0]:
            print(f"WARNING: Atom count mismatch: probe={probe_coords.shape[0]}, ref={ref_coords.shape[0]}")
            return float('inf')

        # Calculate RAW RMSD without alignment
        # Both poses should be in the same coordinate frame from docking
        diff = probe_coords - ref_coords
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
                    'source_file': props['source_file'],
                    'pose_idx': props['pose_idx'],
                    'rmsd': rmsd,
                    'gnina_score': props['gnina_score'],
                    'cnn_score': props['cnn_score'],
                    'cnn_affinity': props['cnn_affinity'],
                    'pass_threshold': rmsd < 2.0  # Standard threshold
                }

                results.append(result)

                print(f"File {props['source_file']} Pose {props['pose_idx']}: RMSD={rmsd:.3f} Å, GNINA={props.get('gnina_score')}")

            except Exception as e:
                print(f"ERROR: Failed to calculate RMSD for pose {props['pose_idx']} in {props['source_file']}: {e}")
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
                'docked_path': str(self.docked_path),
                'total_poses': len(results),
                'passing_poses': sum(1 for r in results if r['pass_threshold']),
                'best_rmsd': results[0]['rmsd'] if results else None,
                'results': results
            }, f, indent=2)

        print(f"Saved JSON report to {json_path}")

        # CSV report
        csv_path = output_path.with_suffix('.csv')
        with open(csv_path, 'w') as f:
            f.write("source_file,pose_idx,rmsd,gnina_score,cnn_score,cnn_affinity,pass_threshold\n")
            for r in results:
                f.write(f"{r['source_file']},{r['pose_idx']},{r['rmsd']:.4f},{r.get('gnina_score')},"
                       f"{r.get('cnn_score')},{r.get('cnn_affinity')},{r['pass_threshold']}\n")

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
                print(f"  {i}. {r['source_file']} Pose {r['pose_idx']}: RMSD={r['rmsd']:.3f} Å, GNINA={r.get('gnina_score')}")
            print("="*60)


def main():
    """CLI entry point."""
    import argparse

    parser = argparse.ArgumentParser(description="Compare docked poses to reference crystal structure")
    parser.add_argument('--receptor', required=True, help='Receptor PDB with reference ligand')
    parser.add_argument('--metadata', required=True, help='Receptor metadata JSON')
    parser.add_argument('--poses', required=False, help='Docked poses SDF file or directory')
    parser.add_argument('--output', required=True, help='Output report path (without extension)')
    
    # Drive arguments
    parser.add_argument('--drive-most-recent', action='store_true', help='Download and analyze the most recent run from Google Drive outputs')
    parser.add_argument('--drive-folder', type=str, help='Download and analyze a specific folder name from Google Drive outputs')
    parser.add_argument('--local-dir', type=str, default='data/outputs', help='Local directory to save/find downloaded Drive folders')

    args = parser.parse_args()
    
    # Handle Drive downloading
    target_poses_path = args.poses
    
    if args.drive_most_recent or args.drive_folder:
        if not drive_auth or not drive_io:
            print("ERROR: Drive modules not available. Cannot fetch from Drive.")
            sys.exit(1)
            
        print("Authenticating to Google Drive...")
        drive_service, _ = drive_auth.setup_drive(verify=False)
        folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]
        
        runs = drive_service.files().list(
            q=f"'{folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
            fields="files(id, name, createdTime)",
            orderBy="createdTime desc"
        ).execute().get("files", [])
        
        if not runs:
            print("ERROR: No runs found in Google Drive outputs folder.")
            sys.exit(1)
            
        target_run = None
        if args.drive_most_recent:
            target_run = runs[0]
        else:
            for run in runs:
                if run['name'] == args.drive_folder:
                    target_run = run
                    break
            if not target_run:
                print(f"ERROR: Run folder '{args.drive_folder}' not found in Google Drive.")
                sys.exit(1)
                
        print(f"Selected Drive run: {target_run['name']} (created: {target_run['createdTime'][:10]})")
        
        # Setup local download path
        project_root = Path(__file__).parent.parent.parent
        local_dir = project_root / args.local_dir / target_run['name']
        
        if not local_dir.exists():
            local_dir.mkdir(parents=True, exist_ok=True)
            print(f"Downloading to {local_dir}...")
            drive_io.download_folder(drive_service, target_run['id'], str(local_dir))
        else:
            print(f"Directory {local_dir} already exists. Skipping download.")
            
        target_poses_path = str(local_dir)

    if not target_poses_path:
        print("ERROR: Must specify either --poses or one of the --drive arguments.")
        parser.print_help()
        sys.exit(1)

    try:
        comparator = PoseComparator(
            receptor_pdb=args.receptor,
            metadata_path=args.metadata,
            docked_path=target_poses_path
        )

        results = comparator.compare_all_poses()
        comparator.generate_report(results, Path(args.output))

        print("Pose comparison complete!")

    except Exception as e:
        print(f"ERROR: Pose comparison failed: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
