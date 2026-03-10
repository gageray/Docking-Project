#!/usr/bin/env python3
"""
Prepare Heavy-Atom Files for Docking
Strips hydrogen atoms from receptor and ligand using OpenBabel and RDKit.
"""

import argparse
import subprocess
import sys
from pathlib import Path
from rdkit import Chem


def strip_receptor_hydrogens(input_pdbqt: Path, output_pdbqt: Path):
    """
    Strip hydrogens from receptor PDBQT using OpenBabel.

    Args:
        input_pdbqt: Input receptor PDBQT with hydrogens
        output_pdbqt: Output receptor PDBQT without hydrogens
    """
    print(f"Stripping hydrogens from receptor: {input_pdbqt}")

    # Use OpenBabel with -d flag (delete hydrogens)
    cmd = [
        "obabel",
        "-ipdbqt", str(input_pdbqt),
        "-opdbqt",
        "-O", str(output_pdbqt),
        "-d"  # Delete hydrogens
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"OpenBabel output: {result.stdout}")

        # Count atoms before/after
        with open(input_pdbqt) as f:
            input_atoms = sum(1 for line in f if line.startswith('ATOM'))

        with open(output_pdbqt) as f:
            output_atoms = sum(1 for line in f if line.startswith('ATOM'))

        print(f"Receptor: {input_atoms} atoms → {output_atoms} heavy atoms ({input_atoms - output_atoms} H removed)")

    except subprocess.CalledProcessError as e:
        print(f"ERROR: OpenBabel failed: {e.stderr}")
        sys.exit(1)


def strip_ligand_hydrogens(input_sdf: Path, output_sdf: Path):
    """
    Strip hydrogens from ligand SDF using RDKit.

    Args:
        input_sdf: Input ligand SDF with hydrogens
        output_sdf: Output ligand SDF without hydrogens
    """
    print(f"\nStripping hydrogens from ligand: {input_sdf}")

    supplier = Chem.SDMolSupplier(str(input_sdf), removeHs=False, sanitize=True)
    writer = Chem.SDWriter(str(output_sdf))

    total_input_atoms = 0
    total_output_atoms = 0

    for idx, mol in enumerate(supplier):
        if mol is None:
            print(f"WARNING: Could not parse molecule {idx}")
            continue

        # Count atoms before
        input_atoms = mol.GetNumAtoms()
        total_input_atoms += input_atoms

        # Remove hydrogens
        mol_heavy = Chem.RemoveHs(mol)

        # Count atoms after
        output_atoms = mol_heavy.GetNumAtoms()
        total_output_atoms += output_atoms

        # Copy properties from original mol
        for prop in mol.GetPropNames():
            mol_heavy.SetProp(prop, mol.GetProp(prop))

        writer.write(mol_heavy)

    writer.close()

    print(f"Ligand: {total_input_atoms} atoms → {total_output_atoms} heavy atoms ({total_input_atoms - total_output_atoms} H removed)")


def main():
    parser = argparse.ArgumentParser(description="Strip hydrogens from receptor and ligand")
    parser.add_argument('--receptor', type=Path, help='Input receptor PDBQT')
    parser.add_argument('--receptor-out', type=Path, help='Output receptor PDBQT (heavy atoms only)')
    parser.add_argument('--ligand', type=Path, help='Input ligand SDF')
    parser.add_argument('--ligand-out', type=Path, help='Output ligand SDF (heavy atoms only)')

    args = parser.parse_args()

    if not (args.receptor or args.ligand):
        print("ERROR: Must specify at least --receptor or --ligand")
        parser.print_help()
        sys.exit(1)

    # Process receptor
    if args.receptor:
        if not args.receptor.exists():
            print(f"ERROR: Receptor not found: {args.receptor}")
            sys.exit(1)

        if not args.receptor_out:
            args.receptor_out = args.receptor.with_name(args.receptor.stem + "_heavy" + args.receptor.suffix)

        args.receptor_out.parent.mkdir(parents=True, exist_ok=True)
        strip_receptor_hydrogens(args.receptor, args.receptor_out)

    # Process ligand
    if args.ligand:
        if not args.ligand.exists():
            print(f"ERROR: Ligand not found: {args.ligand}")
            sys.exit(1)

        if not args.ligand_out:
            args.ligand_out = args.ligand.with_name(args.ligand.stem + "_heavy" + args.ligand.suffix)

        args.ligand_out.parent.mkdir(parents=True, exist_ok=True)
        strip_ligand_hydrogens(args.ligand, args.ligand_out)

    print("\n" + "="*60)
    print("HEAVY-ATOM PREPARATION COMPLETE")
    print("="*60)
    if args.receptor:
        print(f"Receptor: {args.receptor_out}")
    if args.ligand:
        print(f"Ligand: {args.ligand_out}")
    print("="*60)


if __name__ == '__main__':
    main()
