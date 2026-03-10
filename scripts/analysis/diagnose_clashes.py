#!/usr/bin/env python3
"""
Clash Detection for Docked Poses
Calculates steric clashes between receptor and ligand using VdW radii.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Dict, Tuple
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem


# VdW radii (Angstroms) from Bondi
VDW_RADII = {
    'H': 1.20,
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'F': 1.47,
    'P': 1.80,
    'S': 1.80,
    'Cl': 1.75,
    'Br': 1.85,
    'I': 1.98,
}


def load_receptor_from_pdbqt(pdbqt_path: Path) -> Tuple[List[str], np.ndarray, List[Dict]]:
    """
    Load receptor from PDBQT file using direct parsing.

    Returns:
        elements: List of element symbols
        coords: Nx3 array of coordinates
        atom_info: List of atom info dicts
    """
    print(f"Loading receptor from {pdbqt_path}")

    elements = []
    coords = []
    atom_info = []

    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Parse coordinates
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coords.append([x, y, z])

                # Parse element (column 77-78 in PDBQT)
                element_raw = line[77:79].strip() if len(line) > 77 else line[12:14].strip()
                # Clean up element symbol (HD, HS, etc. are all hydrogens in AutoDock atom types)
                if element_raw.startswith('H'):
                    element = 'H'
                elif element_raw.startswith('OA'):
                    element = 'O'
                elif element_raw.startswith('SA'):
                    element = 'S'
                elif len(element_raw) > 1:
                    element = element_raw[0].upper() + element_raw[1:].lower()
                else:
                    element = element_raw.upper()
                elements.append(element)

                # Parse atom info
                residue = line[17:20].strip()
                atom_name = line[12:16].strip()
                res_num = line[22:26].strip()
                chain = line[21].strip()

                atom_info.append({
                    'element': element,
                    'residue': residue,
                    'res_num': res_num,
                    'chain': chain,
                    'atom_name': atom_name,
                    'full_id': f"{residue}{res_num}:{chain}:{atom_name}"
                })

    coords = np.array(coords)
    print(f"Loaded receptor: {len(elements)} atoms")

    return elements, coords, atom_info


def load_ligand(sdf_path: Path) -> Tuple[Chem.Mol, np.ndarray]:
    """Load ligand from SDF file."""
    print(f"Loading ligand from {sdf_path}")

    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mol = supplier[0]

    if mol is None:
        print("ERROR: Could not parse ligand")
        sys.exit(1)

    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])

    print(f"Loaded ligand: {mol.GetNumAtoms()} atoms")

    return mol, coords


def calculate_clashes(
    receptor_elements: List[str],
    receptor_coords: np.ndarray,
    receptor_info: List[Dict],
    ligand_mol: Chem.Mol,
    ligand_coords: np.ndarray,
    clash_threshold: float = 0.6
) -> List[Dict]:
    """
    Calculate all steric clashes between receptor and ligand.

    Args:
        clash_threshold: Distance below VdW sum to consider a clash (Angstroms)

    Returns:
        List of clash dicts with distances and atom info
    """
    print(f"\nCalculating clashes (threshold: {clash_threshold} Å below VdW sum)...")

    clashes = []

    for rec_idx in range(len(receptor_elements)):
        rec_element = receptor_elements[rec_idx]
        rec_coord = receptor_coords[rec_idx]
        rec_vdw = VDW_RADII.get(rec_element, 1.7)
        rec_info = receptor_info[rec_idx] if rec_idx < len(receptor_info) else {'full_id': 'UNKNOWN'}

        for lig_idx in range(ligand_mol.GetNumAtoms()):
            lig_atom = ligand_mol.GetAtomWithIdx(lig_idx)
            lig_element = lig_atom.GetSymbol()
            lig_coord = ligand_coords[lig_idx]
            lig_vdw = VDW_RADII.get(lig_element, 1.7)

            # Calculate distance
            distance = np.linalg.norm(rec_coord - lig_coord)

            # VdW sum
            vdw_sum = rec_vdw + lig_vdw

            # Check for clash
            overlap = vdw_sum - distance

            if overlap > clash_threshold:
                clashes.append({
                    'receptor_atom': rec_info['full_id'],
                    'receptor_element': rec_element,
                    'ligand_atom': f"Lig_{lig_element}{lig_idx}",
                    'ligand_element': lig_element,
                    'distance': float(distance),
                    'vdw_sum': float(vdw_sum),
                    'overlap': float(overlap),
                    'is_h_clash': (rec_element == 'H' or lig_element == 'H')
                })

    # Sort by worst clashes first
    clashes.sort(key=lambda x: x['overlap'], reverse=True)

    return clashes


def analyze_clashes(clashes: List[Dict]) -> Dict:
    """Generate summary statistics for clashes."""
    if not clashes:
        return {
            'total_clashes': 0,
            'h_clashes': 0,
            'heavy_clashes': 0,
            'worst_overlap': 0.0
        }

    h_clashes = [c for c in clashes if c['is_h_clash']]
    heavy_clashes = [c for c in clashes if not c['is_h_clash']]

    return {
        'total_clashes': len(clashes),
        'h_clashes': len(h_clashes),
        'heavy_clashes': len(heavy_clashes),
        'worst_overlap': clashes[0]['overlap'],
        'avg_overlap': np.mean([c['overlap'] for c in clashes]),
        'h_clash_percentage': 100.0 * len(h_clashes) / len(clashes) if clashes else 0.0
    }


def main():
    parser = argparse.ArgumentParser(description="Detect steric clashes between receptor and ligand")
    parser.add_argument('--receptor', type=Path, required=True, help='Receptor PDBQT file')
    parser.add_argument('--ligand', type=Path, required=True, help='Ligand SDF file')
    parser.add_argument('--threshold', type=float, default=0.6, help='Clash threshold (Å below VdW sum)')
    parser.add_argument('--output', type=Path, required=True, help='Output JSON report')
    parser.add_argument('--top-n', type=int, default=20, help='Number of top clashes to report')

    args = parser.parse_args()

    if not args.receptor.exists():
        print(f"ERROR: Receptor not found: {args.receptor}")
        sys.exit(1)

    if not args.ligand.exists():
        print(f"ERROR: Ligand not found: {args.ligand}")
        sys.exit(1)

    # Load molecules
    rec_elements, rec_coords, rec_info = load_receptor_from_pdbqt(args.receptor)
    lig_mol, lig_coords = load_ligand(args.ligand)

    # Calculate clashes
    clashes = calculate_clashes(rec_elements, rec_coords, rec_info, lig_mol, lig_coords, args.threshold)

    # Generate summary
    summary = analyze_clashes(clashes)

    print("\n" + "="*60)
    print("CLASH ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total clashes: {summary['total_clashes']}")
    print(f"H-atom clashes: {summary['h_clashes']} ({summary['h_clash_percentage']:.1f}%)")
    print(f"Heavy-atom clashes: {summary['heavy_clashes']}")
    print(f"Worst overlap: {summary['worst_overlap']:.3f} Å")

    if summary['total_clashes'] > 0:
        print(f"\nTop {min(args.top_n, len(clashes))} worst clashes:")
        for i, clash in enumerate(clashes[:args.top_n], 1):
            h_marker = " [H-CLASH]" if clash['is_h_clash'] else ""
            print(f"  {i}. {clash['receptor_atom']} ({clash['receptor_element']}) <-> "
                  f"{clash['ligand_atom']} ({clash['ligand_element']})")
            print(f"     Distance: {clash['distance']:.3f} Å, VdW sum: {clash['vdw_sum']:.3f} Å, "
                  f"Overlap: {clash['overlap']:.3f} Å{h_marker}")

    print("="*60)

    # Save report
    report = {
        'receptor': str(args.receptor),
        'ligand': str(args.ligand),
        'threshold': args.threshold,
        'summary': summary,
        'top_clashes': clashes[:args.top_n]
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"\nSaved detailed report to {args.output}")


if __name__ == '__main__':
    main()
