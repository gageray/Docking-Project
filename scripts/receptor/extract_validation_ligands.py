"""
Extract validation ligands from ALIGNED receptors.
Must be run AFTER align_model_space.py so ligands are at (0,0,0).
"""
import argparse
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pymol
from pymol import cmd
from utils import setup_logger

logger = setup_logger(__name__)


def extract_ligand_from_aligned(aligned_pdb, resname, output_path):
    """
    Extract ligand from aligned PDB (already centered at origin).

    Args:
        aligned_pdb: Path to aligned PDB file
        resname: Residue name to extract
        output_path: Path to save extracted ligand
    """
    cmd.reinitialize()
    cmd.load(aligned_pdb, "structure")

    # Select ligand
    ligand_sel = f"structure and resn {resname}"
    n_atoms = cmd.count_atoms(ligand_sel)

    if n_atoms == 0:
        logger.info(f"Error: No ligand {resname} found in {aligned_pdb}")
        return False

    logger.info(f"Found {n_atoms} atoms for {resname}")

    # Get centroid to verify it's at origin
    centroid = cmd.centerofmass(ligand_sel)
    logger.info(f"Ligand centroid: ({centroid[0]:.3f}, {centroid[1]:.3f}, {centroid[2]:.3f})")

    # Create clean fragment
    cmd.create("ligand", ligand_sel)
    cmd.alter("ligand", "chain='Z'")
    cmd.alter("ligand", "resi='1'")

    # Save
    cmd.save(output_path, "ligand")
    logger.info(f"Saved to {output_path}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Extract validation ligands from aligned receptors"
    )
    parser.add_argument(
        "--fyp-aligned",
        default="data/receptors/6X3U_aligned.pdb",
        help="Path to aligned 6X3U PDB"
    )
    parser.add_argument(
        "--dzp-aligned",
        default="data/receptors/6X3X_aligned.pdb",
        help="Path to aligned 6X3X PDB"
    )
    parser.add_argument(
        "--output-dir",
        default="data/receptors/validation_ligands",
        help="Output directory"
    )
    args = parser.parse_args()

    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()

    os.makedirs(args.output_dir, exist_ok=True)

    # Extract FYP
    logger.info("\n=== Extracting FYP from 6X3U ===")
    fyp_out = os.path.join(args.output_dir, "FYP_flumazenil.pdb")
    if os.path.exists(args.fyp_aligned):
        success_fyp = extract_ligand_from_aligned(args.fyp_aligned, "FYP", fyp_out)
    else:
        logger.info(f"Error: {args.fyp_aligned} not found. Run align_model_space.py first.")
        success_fyp = False

    # Extract DZP
    logger.info("\n=== Extracting DZP from 6X3X ===")
    dzp_out = os.path.join(args.output_dir, "DZP_diazepam.pdb")
    if os.path.exists(args.dzp_aligned):
        success_dzp = extract_ligand_from_aligned(args.dzp_aligned, "DZP", dzp_out)
    else:
        logger.info(f"Error: {args.dzp_aligned} not found. Run align_model_space.py first.")
        success_dzp = False

    cmd.quit()

    if success_fyp and success_dzp:
        logger.info("\n✓ Successfully extracted validation ligands")
    else:
        logger.info("\n✗ Failed - ensure both receptors are aligned first")
        sys.exit(1)


if __name__ == "__main__":
    main()
