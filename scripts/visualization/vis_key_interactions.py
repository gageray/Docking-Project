#!/usr/bin/env python3
"""
Visualize Key Ligand-Residue Interactions
Compare crystal ligand vs best docked pose, showing only the 2 most important interacting residues for each.
"""

import argparse
import sys
import json
from pathlib import Path

try:
    from pymol import cmd
except ImportError:
    print("ERROR: PyMOL not found. Install with: conda install -c conda-forge pymol-open-source")
    sys.exit(1)


def load_metadata(metadata_path):
    """Load receptor metadata."""
    with open(metadata_path) as f:
        return json.load(f)


def find_closest_residues(ligand_sel, receptor_sel, n=2):
    """
    Find the N closest residues to ligand.
    Returns list of (chain, resi, resn) tuples.
    """
    # Create temporary selection of receptor residues within 5A
    cmd.select("nearby_temp", f"{receptor_sel} within 5 of {ligand_sel}")

    # Get unique residues
    residues = []
    cmd.iterate("nearby_temp and name CA", "residues.append((chain, resi, resn))", space={'residues': residues})

    if not residues:
        return []

    # Calculate minimum distance for each residue
    res_distances = []
    for chain, resi, resn in set(residues):
        res_sel = f"{receptor_sel} and chain {chain} and resi {resi}"
        dist = cmd.distance(f"temp_dist_{resi}", ligand_sel, res_sel, mode=0, quiet=1)
        cmd.delete(f"temp_dist_{resi}")
        if dist > 0:
            res_distances.append((dist, chain, resi, resn))

    # Sort by distance and return top N
    res_distances.sort()
    cmd.delete("nearby_temp")

    return [(chain, resi, resn) for _, chain, resi, resn in res_distances[:n]]


def visualize_key_interactions(receptor_pdb, metadata_path, docked_sdf, output_session=None):
    """
    Create 3-scene visualization:
    1. Crystal ligand + 2 key residues only
    2. Docked ligand + 2 closest residues only
    3. Both ligands superimposed (no protein)
    """
    cmd.reinitialize()

    # Load metadata
    metadata = load_metadata(metadata_path)
    crystal_chain = metadata['target_ligand_chain']
    crystal_resn = metadata['target_ligand_resn']
    key_alpha_res = metadata['key_alpha_res']
    key_gamma_res = metadata['key_gamma_res']

    print("="*80)
    print("LOADING STRUCTURES")
    print("="*80)
    print(f"Crystal ligand: {crystal_resn} chain {crystal_chain}")
    print(f"Key residues from metadata: alpha-{key_alpha_res}, gamma-{key_gamma_res}")

    # Load receptor with crystal ligand
    cmd.load(str(receptor_pdb), "full_structure")

    # Load docked pose
    cmd.load(str(docked_sdf), "docked_ligand")

    # Extract crystal ligand as separate object
    crystal_sel = f"full_structure and chain {crystal_chain} and resn {crystal_resn}"
    cmd.create("crystal_ligand", crystal_sel)

    # Extract just the protein (no ligand)
    cmd.create("receptor", f"full_structure and not (chain {crystal_chain} and resn {crystal_resn})")

    # Delete full structure
    cmd.delete("full_structure")

    print("\n" + "="*80)
    print("IDENTIFYING KEY RESIDUES")
    print("="*80)

    # Find alpha and gamma chains (a1_bzd and y2_bzd from metadata)
    # Based on metadata: D=a1_bzd, E=y2_bzd
    alpha_chain = "D"  # a1_bzd
    gamma_chain = "E"  # y2_bzd

    # Crystal ligand key residues (from metadata)
    crystal_key_res = [
        (alpha_chain, key_alpha_res, ""),
        (gamma_chain, key_gamma_res, "")
    ]

    # Get actual residue names
    crystal_key_residues = []
    for chain, resi, _ in crystal_key_res:
        res_sel = f"receptor and chain {chain} and resi {resi} and name CA"
        if cmd.count_atoms(res_sel) > 0:
            resn = cmd.get_model(res_sel).atom[0].resn
            crystal_key_residues.append((chain, resi, resn))
            print(f"Crystal key residue: {resn}{resi} chain {chain}")

    # Use THE SAME key residues for both (from metadata, not auto-detected)
    # These are the known important binding residues
    print(f"\nUsing SAME key residues for both crystal and docked:")
    for chain, resi, resn in crystal_key_residues:
        print(f"  {resn}{resi} chain {chain}")

    # Extract key residues as ONE object (same residues used for both scenes)
    key_res_sel = " or ".join([f"(chain {chain} and resi {resi})" for chain, resi, _ in crystal_key_residues])
    cmd.create("key_residues", f"receptor and ({key_res_sel})")
    print(f"\nCreated key_residues object: {cmd.count_atoms('key_residues')} atoms")
    print(f"Crystal ligand atoms: {cmd.count_atoms('crystal_ligand')} atoms")
    print(f"Docked ligand atoms: {cmd.count_atoms('docked_ligand')} atoms")

    # DELETE the entire receptor - we only want the extracted residues, no ghost protein
    cmd.delete("receptor")

    # Verify objects exist
    print(f"Objects in session: {cmd.get_object_list()}")

    print("\n" + "="*80)
    print("CREATING VISUALIZATION SCENES")
    print("="*80)

    # ========== SCENE 1: CRYSTAL LIGAND + KEY RESIDUES ONLY ==========
    # Hide everything first
    cmd.hide("everything", "all")

    # Show only crystal ligand and key residues
    cmd.show("sticks", "crystal_ligand")
    cmd.show("sticks", "key_residues")
    cmd.util.cnc("key_residues")  # Color by element with carbons colored
    cmd.color("cyan", "key_residues and elem C")  # Override carbons to cyan

    # Style crystal ligand
    cmd.color("green", "crystal_ligand and elem C")
    cmd.color("red", "crystal_ligand and elem O")
    cmd.color("blue", "crystal_ligand and elem N")
    cmd.color("cyan", "crystal_ligand and elem F")
    cmd.set("stick_radius", 0.2, "crystal_ligand")

    # Style key residues - make sure all atoms are colored
    cmd.set("stick_radius", 0.2, "key_residues")

    # Add labels
    for chain, resi, resn in crystal_key_residues:
        cmd.label(f"key_residues and chain {chain} and resi {resi} and name CA", f'"{resn}{resi}"')

    cmd.center("crystal_ligand")
    cmd.zoom("crystal_ligand or key_residues", buffer=5)
    cmd.bg_color("white")

    cmd.scene("crystal_only", "store", message="Crystal ligand + key residues (PHE101, PHE77)")
    print("Scene 1: crystal_only (Crystal ligand GREEN + key residues CYAN)")

    # ========== SCENE 2: DOCKED LIGAND + SAME KEY RESIDUES ==========
    # Hide everything first
    cmd.hide("everything", "all")

    # Show only docked ligand and SAME key residues
    cmd.show("sticks", "docked_ligand")
    cmd.show("sticks", "key_residues")
    cmd.util.cnc("key_residues")  # Color by element
    cmd.color("orange", "key_residues and elem C")  # Override carbons to orange

    # Style docked ligand
    cmd.color("magenta", "docked_ligand and elem C")
    cmd.color("red", "docked_ligand and elem O")
    cmd.color("blue", "docked_ligand and elem N")
    cmd.color("cyan", "docked_ligand and elem F")
    cmd.set("stick_radius", 0.2, "docked_ligand")

    # Style key residues - make sure all atoms are colored
    cmd.set("stick_radius", 0.2, "key_residues")

    # Add labels
    for chain, resi, resn in crystal_key_residues:
        cmd.label(f"key_residues and chain {chain} and resi {resi} and name CA", f'"{resn}{resi}"')

    cmd.center("docked_ligand")
    cmd.zoom("docked_ligand or key_residues", buffer=5)

    cmd.scene("docked_only", "store", message="Docked ligand + closest residues")
    print("Scene 2: docked_only (Docked ligand MAGENTA + closest residues ORANGE)")

    # ========== SCENE 3: BOTH LIGANDS SUPERIMPOSED (NO PROTEIN) ==========
    # Hide everything first
    cmd.hide("everything", "all")

    # Show only both ligands
    cmd.show("sticks", "crystal_ligand")
    cmd.show("sticks", "docked_ligand")

    # Style crystal ligand
    cmd.color("green", "crystal_ligand and elem C")
    cmd.color("red", "crystal_ligand and elem O")
    cmd.color("blue", "crystal_ligand and elem N")
    cmd.color("cyan", "crystal_ligand and elem F")
    cmd.set("stick_radius", 0.18, "crystal_ligand")

    # Style docked ligand
    cmd.color("magenta", "docked_ligand and elem C")
    cmd.color("red", "docked_ligand and elem O")
    cmd.color("blue", "docked_ligand and elem N")
    cmd.color("cyan", "docked_ligand and elem F")
    cmd.set("stick_radius", 0.18, "docked_ligand")

    # Add labels
    cmd.label("crystal_ligand and name C1", '"Crystal"')
    cmd.label("docked_ligand and name C1", '"Docked"')
    cmd.set("label_color", "green", "crystal_ligand")
    cmd.set("label_color", "magenta", "docked_ligand")

    cmd.center("crystal_ligand")
    cmd.zoom("crystal_ligand or docked_ligand", buffer=5)

    cmd.scene("superimposed", "store", message="Both ligands superimposed (no protein)")
    print("Scene 3: superimposed (Crystal GREEN vs Docked MAGENTA)")

    # Calculate RMSD
    print("\n" + "="*80)
    print("RMSD CALCULATION")
    print("="*80)
    try:
        # Try aligned RMSD
        rmsd_aligned = cmd.align("docked_ligand and name C*", "crystal_ligand and name C*")[0]
        print(f"Aligned RMSD (after superposition): {rmsd_aligned:.3f} Å")

        # Calculate unaligned RMSD (reset first)
        cmd.load(str(docked_sdf), "docked_temp")
        rmsd_unaligned = cmd.rms_cur("docked_temp and name C*", "crystal_ligand and name C*")
        print(f"Unaligned RMSD (absolute coordinates): {rmsd_unaligned:.3f} Å")
        cmd.delete("docked_temp")
    except Exception as e:
        print(f"Could not calculate RMSD: {e}")

    # Save session
    if output_session:
        cmd.save(str(output_session))
        print(f"\n✓ Saved PyMOL session to: {output_session}")

    print("\n" + "="*80)
    print("VISUALIZATION COMPLETE")
    print("="*80)
    print("To switch between scenes in PyMOL:")
    print("  scene crystal_only    - Crystal ligand + key residues")
    print("  scene docked_only     - Docked ligand + closest residues")
    print("  scene superimposed    - Both ligands overlaid")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize key ligand-residue interactions"
    )
    parser.add_argument(
        "--receptor",
        type=Path,
        required=True,
        help="Receptor PDB with crystal ligand"
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Receptor metadata JSON"
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

    args = parser.parse_args()

    if not args.receptor.exists():
        print(f"ERROR: Receptor not found: {args.receptor}")
        sys.exit(1)

    if not args.metadata.exists():
        print(f"ERROR: Metadata not found: {args.metadata}")
        sys.exit(1)

    if not args.docked.exists():
        print(f"ERROR: Docked pose not found: {args.docked}")
        sys.exit(1)

    visualize_key_interactions(
        args.receptor,
        args.metadata,
        args.docked,
        args.output
    )


if __name__ == "__main__":
    main()
