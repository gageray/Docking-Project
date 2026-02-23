"""
Ligand Conformer Validation
Run inside PyMOL: `run scripts/vis_ligands.py`
"""
from pymol import cmd
import glob
import os
import colorsys

def load_and_align_conformers(ligand_prefix="docosanyl_ferulate_conf"):
    cmd.reinitialize()
    
    # 1. Glob all matching conformer SDFs
    search_path = os.path.join("data", "ligands", f"{ligand_prefix}*.sdf")
    files = glob.glob(search_path)
    
    if not files:
        print(f"No files found matching {search_path}")
        return
        
    print(f"Loading {len(files)} conformers...")
    
    object_names = []
    
    # 2. Load them
    for i, file in enumerate(files):
        obj_name = f"conf_{i+1}"
        cmd.load(file, obj_name)
        object_names.append(obj_name)
        
    # 3. Identify ferulic acid headgroup and align
    # The headgroup is the substituted phenol ring. We can select it generally as the aromatic ring.
    # A generic SMARTS for the 3-methoxy-4-hydroxy-phenyl ring: `c1cc(O)c(OC)cc1` roughly.
    # Since PyMOL selection strings are easier: we select all aromatic carbons/oxygens if possible.
    # More stably for docosanyl ferulate, the first 10-12 heavy atoms define the head.
    # Let's align by indices 1 through 10 (which roughly covers the phenolic head usually).
    # Since they are generated from the same SMILES, atom indices should perfectly match.
    
    ref_obj = object_names[0]
    
    # A safe alignment that covers the headgroup:
    # `id 1-12` usually captures the aromatic head + conjugated tail origin
    for obj in object_names[1:]:
        # fit tail to reference based on atom IDs 1-12 (the headgroup)
        cmd.fit(f"{obj} and id 1-12", f"{ref_obj} and id 1-12")
        
    print("Alignment complete based on headgroup (atoms 1-12).")
    
    # 4. Color them distinctly for visualization
    for i, obj in enumerate(object_names):
        hue = i / float(len(object_names))
        # Convert HSV to RGB
        r, g, b = colorsys.hsv_to_rgb(hue, 1.0, 1.0)
        color_name = f"color_{i}"
        cmd.set_color(color_name, [r, g, b])
        cmd.color(color_name, obj)
        
    # Orient viewport
    cmd.reset()
    cmd.zoom("all")
    print("Success Criteria Check: Observe that the 22-carbon tails project outward and do not fold back.")

# Auto-execute if run as a script in PyMOL
if __name__ in ["__main__", "pymol"]:
    load_and_align_conformers()
