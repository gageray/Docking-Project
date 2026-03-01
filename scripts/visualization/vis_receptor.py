"""
Receptor Visualization Script
Generates publication-quality images of aligned receptor structures.
"""
import argparse
import json
import pymol
from pymol import cmd


def setup_scene(obj, meta):
    """Apply consistent color scheme and display settings."""
    cmd.bg_color("white")
    cmd.hide("all")
    cmd.show("cartoon", obj)

    # Color scheme: a1/y2_bzd colored, others grey
    chains = meta.get("chains", {})
    for chain_id, chain_desc in chains.items():
        if "_bzd" in chain_desc:
            # Color BZD pocket chains
            if "a1" in chain_desc or "alpha1" in chain_desc:
                cmd.color("cyan", f"{obj} and chain {chain_id}")
                # Highlight key alpha residue
                key_res = meta.get("key_alpha_res", "101")
                cmd.color("purple", f"{obj} and chain {chain_id} and resi {key_res}")
                cmd.show("sticks", f"{obj} and chain {chain_id} and resi {key_res}")
            elif "y2" in chain_desc or "gamma2" in chain_desc:
                cmd.color("yellow", f"{obj} and chain {chain_id}")
                # Highlight key gamma residue
                key_res = meta.get("key_gamma_res", "77")
                cmd.color("purple", f"{obj} and chain {chain_id} and resi {key_res}")
                cmd.show("sticks", f"{obj} and chain {chain_id} and resi {key_res}")
            else:
                cmd.color("magenta", f"{obj} and chain {chain_id}")  # Fallback for other bzd
        else:
            # Grey for non-BZD chains
            cmd.color("grey70", f"{obj} and chain {chain_id}")

    # Color ligand if present (bright green carbons, element-specific for others)
    ligand_resn = meta.get("target_ligand_resn")
    ligand_chain = meta.get("target_ligand_chain")
    if ligand_resn and ligand_chain:
        lig_sel = f"{obj} and resn {ligand_resn} and chain {ligand_chain}"
        if cmd.count_atoms(lig_sel) > 0:
            cmd.show("sticks", lig_sel)
            cmd.color("green", f"{lig_sel} and elem C")
            cmd.color("red", f"{lig_sel} and elem O")
            cmd.color("red", f"{lig_sel} and elem N")
            cmd.color("slate", f"{lig_sel} and elem F")
            cmd.color("grey90", f"{lig_sel} and elem H")


def render_side_view(obj, output_path, width=1920, height=1080):
    """Full frame side view of receptor."""
    cmd.zoom(obj, complete=1)
    cmd.turn('x', 270)
    cmd.turn('y', 270)
    cmd.png(output_path, width=width, height=height, dpi=300, ray=1)

def render_pocket_view(obj, meta, output_path, width=1920, height=1080):
    """Zoomed view of BZD binding pocket."""
    # Apply user-requested transparency for pocket view
    # Note: PyMOL transparency is 0.0 (opaque) to 1.0 (invisible)
    # User wanted "alpha 0 for grey", "alpha 60 for BZD", "alpha 100 for key/lig"
    # Mapping alpha logic: alpha 0% = transparency 1.0, alpha 60% = transparency 0.4, alpha 100% = transparency 0.0
    
    cmd.set("cartoon_transparency", 1.0, f"{obj}") # Make all chains invisible by default
    
    # Identify bzd chains and make them alpha 60% (transparency 0.4)
    chains = meta.get("chains", {})
    bzd_chains = [c for c, desc in chains.items() if "_bzd" in desc]
    if bzd_chains:
        bzd_sel = " or ".join([f"chain {c}" for c in bzd_chains])
        cmd.set("cartoon_transparency", 0.4, f"{obj} and ({bzd_sel})")
    
    # Leave stick transparency at 0.0 (alpha 100%) which applies to the ligand and key residues
    cmd.set("stick_transparency", 0.0, obj)

    ligand_resn = meta.get("target_ligand_resn")
    ligand_chain = meta.get("target_ligand_chain")

    if ligand_resn and ligand_chain:
        lig_sel = f"{obj} and resn {ligand_resn} and chain {ligand_chain}"
        if cmd.count_atoms(lig_sel) > 0:
            # Zoom to pocket area (ligand + 8Å around it)
            cmd.zoom(f"{lig_sel} or ({obj} within 8 of {lig_sel})", buffer=3)
            cmd.png(output_path, width=width, height=height, dpi=300, ray=1)
            return

    # No ligand (or missing), zoom to the key residues instead of a generic center
    key_alpha = meta.get("key_alpha_res", "101")
    key_gamma = meta.get("key_gamma_res", "77")
    
    alpha_chain = None
    gamma_chain = None
    for chain_id, desc in chains.items():
        if "_bzd" in desc:
            if "a1" in desc or "alpha1" in desc:
                alpha_chain = chain_id
            elif "y2" in desc or "gamma2" in desc:
                gamma_chain = chain_id
                
    sel_parts = []
    if alpha_chain: sel_parts.append(f"(chain {alpha_chain} and resi {key_alpha})")
    if gamma_chain: sel_parts.append(f"(chain {gamma_chain} and resi {key_gamma})")
    
    if sel_parts:
        key_sel = f"{obj} and ({' or '.join(sel_parts)})"
        cmd.zoom(f"{key_sel} or ({obj} within 10 of {key_sel})", buffer=3)
    else:
        # Fallback
        cmd.pseudoatom("pocket_center", pos=[0.0, 0.0, 0.0])
        cmd.zoom(f"{obj} within 10 of pocket_center", buffer=1)

    cmd.png(output_path, width=width, height=height, dpi=300, ray=1)


def main():
    parser = argparse.ArgumentParser(description="Visualize aligned receptor structures")
    parser.add_argument("--pdb", required=True, help="Aligned receptor PDB file")
    parser.add_argument("--json", required=True, help="Receptor metadata JSON")
    parser.add_argument("--out_prefix", required=True, help="Output prefix for PNG files")
    parser.add_argument("--width", type=int, default=1920, help="Image width (default: 1920)")
    parser.add_argument("--height", type=int, default=1080, help="Image height (default: 1080)")
    args = parser.parse_args()

    # Initialize PyMOL
    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()
    cmd.reinitialize()

    # Load metadata
    with open(args.json, 'r') as f:
        meta = json.load(f)

    # Load structure
    cmd.load(args.pdb, "receptor")

    # Setup visualization
    setup_scene("receptor", meta)

    # Render images
    render_side_view("receptor", f"{args.out_prefix}_side.png", args.width, args.height)
    render_pocket_view("receptor", meta, f"{args.out_prefix}_pocket.png", args.width, args.height)

    cmd.quit()
    print(f"Generated: {args.out_prefix}_side.png")
    print(f"Generated: {args.out_prefix}_pocket.png")


if __name__ == "__main__":
    main()
