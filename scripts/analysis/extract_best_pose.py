#!/usr/bin/env python3
"""
Extract Best Pose from Docking Results
Extracts a specific pose from GNINA output SDF for further analysis.
"""

import argparse
import sys
from pathlib import Path
from rdkit import Chem


def extract_pose(input_sdf: Path, pose_index: int, output_sdf: Path):
    """
    Extract a specific pose from multi-pose SDF file.

    Args:
        input_sdf: Path to GNINA output SDF with multiple poses
        pose_index: 0-based index of pose to extract
        output_sdf: Path for output SDF file
    """
    print(f"Loading poses from {input_sdf}")

    supplier = Chem.SDMolSupplier(str(input_sdf), removeHs=False, sanitize=True)

    # Count total poses
    mols = [mol for mol in supplier if mol is not None]
    total_poses = len(mols)

    print(f"Found {total_poses} valid poses in SDF")

    if pose_index >= total_poses or pose_index < 0:
        print(f"ERROR: Pose index {pose_index} out of range [0-{total_poses-1}]")
        sys.exit(1)

    # Extract target pose
    target_mol = mols[pose_index]

    # Print properties
    print(f"\nExtracting pose {pose_index}:")
    props = target_mol.GetPropsAsDict()
    for key, value in props.items():
        print(f"  {key}: {value}")

    # Write to output
    writer = Chem.SDWriter(str(output_sdf))
    writer.write(target_mol)
    writer.close()

    print(f"\nSaved pose {pose_index} to {output_sdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract specific pose from GNINA docking output"
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Input SDF file with multiple poses"
    )
    parser.add_argument(
        "--pose-index",
        type=int,
        required=True,
        help="0-based index of pose to extract"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output SDF file for extracted pose"
    )

    args = parser.parse_args()

    if not args.input.exists():
        print(f"ERROR: Input file not found: {args.input}")
        sys.exit(1)

    # Create output directory if needed
    args.output.parent.mkdir(parents=True, exist_ok=True)

    extract_pose(args.input, args.pose_index, args.output)


if __name__ == "__main__":
    main()
