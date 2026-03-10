#!/usr/bin/env python3
"""
Visualize Crystal Ligand vs Best Docked Pose
Loads receptor with crystal ligand (chain Z) and overlays best docked pose for comparison.
"""

import argparse
import sys
from pathlib import Path

try:
    from pymol import cmd
except ImportError:
    print("ERROR: PyMOL not found. Install with: conda install -c conda-forge pymol-open-source")
    sys.exit(1)


def visualize_comparison(receptor_pdb: Path, crystal_chain: str, crystal_resn: str,
                         docked_sdf: Path, output_session: Path = None):
    """
    Create PyMOL visualization comparing crystal ligand to docked pose.

    Args:
        receptor_pdb: PDB file with crystal structure (must contain ligand in chain Z)
        crystal_chain: Chain ID of crystal ligand
        crystal_resn: Residue name of crystal ligand
        docked_sdf: SDF file of best docked pose
        output_session: Optional path to save PyMOL session
    """
    cmd.reinitialize()

    # Load receptor with crystal ligand
    cmd.load(str(receptor_pdb), "receptor")

    # Load docked pose
    cmd.load(str(docked_sdf), "docked_pose")

    # Extract crystal ligand
    crystal_sel = f"receptor and chain {crystal_chain} and resn {crystal_resn}"
    cmd.create("crystal_ligand", crystal_sel)

    # Remove ligand from receptor display
    cmd.remove(f"receptor and chain {crystal_chain}")

    # Style receptor
    cmd.hide("everything", "receptor")
    cmd.show("cartoon", "receptor")
    cmd.color("gray70", "receptor")

    # Show binding pocket residues (within 5A of crystal ligand)
    cmd.select("pocket", f"receptor within 5 of crystal_ligand")
    cmd.show("sticks", "pocket")
    cmd.color("wheat", "pocket and elem C")

    # Style crystal ligand
    cmd.hide("everything", "crystal_ligand")
    cmd.show("sticks", "crystal_ligand")
    cmd.color("green", "crystal_ligand and elem C")
    cmd.color("red", "crystal_ligand and elem O")
    cmd.color("blue", "crystal_ligand and elem N")
    cmd.color("cyan", "crystal_ligand and elem F")
    cmd.set("stick_radius", 0.15, "crystal_ligand")

    # Style docked pose
    cmd.hide("everything", "docked_pose")
    cmd.show("sticks", "docked_pose")
    cmd.color("magenta", "docked_pose and elem C")
    cmd.color("red", "docked_pose and elem O")
    cmd.color("blue", "docked_pose and elem N")
    cmd.color("cyan", "docked_pose and elem F")
    cmd.set("stick_radius", 0.15, "docked_pose")

    # Add labels
    cmd.label("crystal_ligand and name C1", '"Crystal"')
    cmd.label("docked_pose and name C1", '"Docked"')
    cmd.set("label_color", "green", "crystal_ligand")
    cmd.set("label_color", "magenta", "docked_pose")

    # Calculate RMSD
    try:
        rmsd = cmd.rms_cur("docked_pose and name C*", "crystal_ligand and name C*")
        print(f"\nHeavy atom RMSD: {rmsd:.3f} Å")
    except:
        print("Could not calculate RMSD (atom count mismatch)")

    # Measure distances between closest atoms
    cmd.distance("clash_check", "docked_pose", "pocket", cutoff=2.5, mode=2)
    cmd.hide("labels", "clash_check")
    cmd.color("red", "clash_check")

    # Center on ligands
    cmd.center("crystal_ligand")
    cmd.zoom("crystal_ligand or docked_pose", buffer=5)

    # Set background
    cmd.bg_color("white")

    # Show origin sphere to verify alignment
    cmd.pseudoatom("origin", pos=[0, 0, 0])
    cmd.show("sphere", "origin")
    cmd.color("yellow", "origin")
    cmd.set("sphere_scale", 1.0, "origin")

    print("\n" + "="*60)
    print("VISUALIZATION LOADED")
    print("="*60)
    print("Objects:")
    print("  - receptor (gray cartoon + pocket residues)")
    print(f"  - crystal_ligand (GREEN, chain {crystal_chain}, {crystal_resn})")
    print("  - docked_pose (MAGENTA)")
    print("  - origin (YELLOW sphere at 0,0,0)")
    print("  - clash_check (red distances < 2.5 Å)")
    print("="*60)

    # Save session if requested
    if output_session:
        cmd.save(str(output_session))
        print(f"\nSaved PyMOL session to: {output_session}")

    # Print clashing residues
    print("\nPocket residues within 5 Å of crystal ligand:")
    cmd.iterate("pocket and name CA", "print(f'{resn}{resi}:{chain}')")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize crystal ligand vs docked pose in PyMOL"
    )
    parser.add_argument(
        "--receptor",
        type=Path,
        required=True,
        help="Receptor PDB with crystal ligand"
    )
    parser.add_argument(
        "--crystal-chain",
        type=str,
        default="Z",
        help="Chain ID of crystal ligand (default: Z)"
    )
    parser.add_argument(
        "--crystal-resn",
        type=str,
        default="FYP",
        help="Residue name of crystal ligand (default: FYP)"
    )
    parser.add_argument(
        "--docked",
        type=Path,
        required=True,
        help="Best docked pose SDF file"
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output PyMOL session file (.pse)"
    )
    parser.add_argument(
        "--interactive",
        action="store_true",
        help="Keep PyMOL window open (requires PyMOL GUI)"
    )

    args = parser.parse_args()

    if not args.receptor.exists():
        print(f"ERROR: Receptor not found: {args.receptor}")
        sys.exit(1)

    if not args.docked.exists():
        print(f"ERROR: Docked pose not found: {args.docked}")
        sys.exit(1)

    visualize_comparison(
        args.receptor,
        args.crystal_chain,
        args.crystal_resn,
        args.docked,
        args.output
    )

    if args.interactive:
        print("\nLaunching PyMOL GUI...")
        cmd.do("_ cmd.window('show')")
    else:
        print("\nTo view interactively, use --interactive flag")


if __name__ == "__main__":
    main()
