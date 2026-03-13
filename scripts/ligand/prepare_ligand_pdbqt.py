#!/usr/bin/env python3
"""
Convert ligands from MOL2 format to PDBQT using Meeko.
Preserves aromatic bond types for accurate GNINA scoring.

Usage:
    python prepare_ligand_pdbqt.py --input-dir data/ligands --output-dir data/ligands/pdbqt
    python prepare_ligand_pdbqt.py --input-file ligand.mol2 --output-file ligand.pdbqt
"""

import argparse
import subprocess
import zipfile
import shutil
from pathlib import Path
import sys


def convert_mol2_to_pdbqt(mol2_path: Path, pdbqt_path: Path) -> bool:
    """
    Convert a single MOL2 file to PDBQT using Meeko.

    Args:
        mol2_path: Path to input MOL2 file
        pdbqt_path: Path to output PDBQT file

    Returns:
        True if conversion succeeded, False otherwise
    """
    # Meeko's mk_prepare_ligand.py command
    # -i: input file
    # -o: output file (without .pdbqt extension)
    out_basename = str(pdbqt_path.with_suffix(''))

    cmd = [
        "mk_prepare_ligand.py",
        "-i", str(mol2_path),
        "-o", out_basename
    ]

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        print(f"  ✓ Converted {mol2_path.name} → {pdbqt_path.name}")
        return True

    except subprocess.CalledProcessError as e:
        print(f"  ✗ Failed to convert {mol2_path.name}")
        if e.stderr:
            print(f"    Error: {e.stderr}")
        return False

    except FileNotFoundError:
        print("  ✗ ERROR: mk_prepare_ligand.py not found in PATH")
        print("    Install Meeko: conda install -c conda-forge meeko")
        sys.exit(1)


def process_zip_archive(zip_path: Path, output_dir: Path) -> int:
    """
    Extract MOL2 files from zip archive and convert to PDBQT.

    Args:
        zip_path: Path to zip file containing MOL2 files
        output_dir: Directory to save PDBQT files

    Returns:
        Number of files successfully converted
    """
    temp_dir = output_dir / f"_temp_{zip_path.stem}"
    temp_dir.mkdir(parents=True, exist_ok=True)

    converted_count = 0

    try:
        # Extract zip to temporary directory
        with zipfile.ZipFile(zip_path, 'r') as zf:
            mol2_files = [name for name in zf.namelist() if name.endswith('.mol2')]

            if not mol2_files:
                print(f"  ! No MOL2 files found in {zip_path.name}")
                return 0

            print(f"[*] Processing {zip_path.name} ({len(mol2_files)} conformers)")

            for mol2_name in mol2_files:
                # Extract to temp directory
                zf.extract(mol2_name, temp_dir)
                mol2_path = temp_dir / mol2_name

                # Convert to PDBQT
                pdbqt_name = mol2_name.replace('.mol2', '.pdbqt')
                pdbqt_path = output_dir / pdbqt_name

                if convert_mol2_to_pdbqt(mol2_path, pdbqt_path):
                    converted_count += 1

    finally:
        # Clean up temporary directory
        if temp_dir.exists():
            shutil.rmtree(temp_dir)

    return converted_count


def process_directory(input_dir: Path, output_dir: Path, recursive: bool = False) -> None:
    """
    Process all MOL2 files and zip archives in a directory.

    Args:
        input_dir: Directory containing MOL2 files and/or zip archives
        output_dir: Directory to save PDBQT files
        recursive: Search subdirectories
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    total_converted = 0

    # Find MOL2 files
    pattern = "**/*.mol2" if recursive else "*.mol2"
    mol2_files = list(input_dir.glob(pattern))

    if mol2_files:
        print(f"[*] Found {len(mol2_files)} MOL2 files")
        for mol2_path in mol2_files:
            pdbqt_name = mol2_path.stem + '.pdbqt'
            pdbqt_path = output_dir / pdbqt_name

            if convert_mol2_to_pdbqt(mol2_path, pdbqt_path):
                total_converted += 1

    # Find zip archives
    pattern = "**/*.zip" if recursive else "*.zip"
    zip_files = list(input_dir.glob(pattern))

    if zip_files:
        print(f"\n[*] Found {len(zip_files)} zip archives")
        for zip_path in zip_files:
            count = process_zip_archive(zip_path, output_dir)
            total_converted += count

    print(f"\n[+] Conversion complete: {total_converted} ligands converted to PDBQT")


def main():
    parser = argparse.ArgumentParser(
        description="Convert MOL2 ligands to PDBQT using Meeko",
        epilog="Example: python prepare_ligand_pdbqt.py --input-dir data/ligands --output-dir data/ligands/pdbqt"
    )

    parser.add_argument(
        '--input-dir',
        type=Path,
        help='Directory containing MOL2 files and/or zip archives'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        help='Directory to save PDBQT files'
    )
    parser.add_argument(
        '--input-file',
        type=Path,
        help='Single MOL2 file to convert'
    )
    parser.add_argument(
        '--output-file',
        type=Path,
        help='Output PDBQT file path (for single file conversion)'
    )
    parser.add_argument(
        '--recursive',
        action='store_true',
        help='Search subdirectories for MOL2 files'
    )

    args = parser.parse_args()

    # Validate arguments
    if args.input_file:
        # Single file mode
        if not args.output_file:
            print("ERROR: --output-file required when using --input-file")
            sys.exit(1)

        if not args.input_file.exists():
            print(f"ERROR: Input file not found: {args.input_file}")
            sys.exit(1)

        args.output_file.parent.mkdir(parents=True, exist_ok=True)
        convert_mol2_to_pdbqt(args.input_file, args.output_file)

    elif args.input_dir:
        # Directory mode
        if not args.output_dir:
            print("ERROR: --output-dir required when using --input-dir")
            sys.exit(1)

        if not args.input_dir.exists():
            print(f"ERROR: Input directory not found: {args.input_dir}")
            sys.exit(1)

        process_directory(args.input_dir, args.output_dir, args.recursive)

    else:
        print("ERROR: Either --input-dir or --input-file must be specified")
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
