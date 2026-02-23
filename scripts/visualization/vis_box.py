"""
Receptor Alignment and Box Centroid
Run inside PyMOL: `run scripts/vis_box.py`
"""
from pymol import cmd
import os
import json

def visualize_boxes():
    cmd.reinitialize()
    
    config_path = "config.json"
    if not os.path.exists(config_path):
        print("config.json not found in root")
        return
        
    try:
        with open(config_path, "r") as f:
            config = json.load(f)
        
        box_cfg = config.get("box_definition", {})
        reference_cif = box_cfg.get("reference_box", {}).get("raw_cif")
        targets = box_cfg.get("targets", [])
        box_params_file = box_cfg.get("output_file", "box_params.json")
        
        if not os.path.exists(box_params_file):
            print(f"Box params file {box_params_file} not found. Run scripts/box_definition.py first!")
            return
            
        with open(box_params_file, "r") as f:
            box_params = json.load(f)
            
    except Exception as e:
        print(f"Error parsing JSON info: {e}")
        return
        
    # 1. Load Reference CIF
    if reference_cif and os.path.exists(reference_cif):
        cmd.load(reference_cif, "reference_crystal")
        cmd.hide("everything", "reference_crystal")
        cmd.show("cartoon", "reference_crystal")
        # Extract and highlight FYP or DZP
        cmd.select("native_ligand", "reference_crystal and (resn DZP or resn FYP)")
        cmd.show("sticks", "native_ligand")
        cmd.color("green", "native_ligand")
        print(f"Loaded {reference_cif}")
    
    # 2. Iterate targets and load PDBs + add Pseudoatoms
    # Keep track of objects
    target_objs = []
    
    for t_data in box_params.get("targets", []):
        t_name = t_data.get("name")
        t_centroid = t_data.get("centroid")
        
        # We need the prepped pdb path from config
        t_pdb = next((t["prepped_pdb"] for t in targets if t["name"] == t_name), None)
        
        if t_pdb and os.path.exists(t_pdb):
            cmd.load(t_pdb, t_name)
            cmd.hide("everything", t_name)
            cmd.show("cartoon", t_name)
            cmd.color("cyan", t_name)
            target_objs.append(t_name)
        else:
            print(f"Could not load PDB for {t_name}")
            
        if t_centroid:
            # Create Pseudoatom
            x, y, z = t_centroid
            p_name = f"box_center_{t_name}"
            cmd.pseudoatom(p_name, pos=[x, y, z])
            cmd.show("spheres", p_name)
            cmd.color("red", p_name)
            cmd.set("sphere_scale", 1.5, p_name)
            
            print(f"Placed {p_name} at {x:.2f}, {y:.2f}, {z:.2f}")

    # Ensure reference centroid is also graphed
    ref_param = box_params.get("reference", {})
    if ref_param.get("centroid"):
        x, y, z = ref_param.get("centroid")
        cmd.pseudoatom("box_center_ref", pos=[x, y, z])
        cmd.show("spheres", "box_center_ref")
        cmd.color("yellow", "box_center_ref")
        cmd.set("sphere_scale", 1.5, "box_center_ref")
    
    # 3. Align everything structurally to the target space for viewing (or reference)
    # The box mapping script calculated target centroids without moving the receptor!
    # Therefore, we do NOT run `cmd.align` because the bounding box points are calculated 
    # to sit precisely in the native target coordinate frame! Aligning them here would decouple 
    # the target from its generated pseudoatom.
    # 
    # For user sanity, we instead align the REFERENCE structure TO the targets to see if the
    # native ligand overlaps the pseudoatoms. We align ref -> target 1
    
    if target_objs and "reference_crystal" in cmd.get_names("objects"):
        print("Aligning reference crystal to first target to verify box overlap...")
        cmd.align("reference_crystal", target_objs[0])
    
    cmd.reset()
    cmd.zoom("box_center_*", buffer=15.0)
    
    print("Success Criteria Check: Ensure the red pseudoatoms sit precisely in the target binding pockets.")

if __name__ in ["__main__", "pymol"]:
    visualize_boxes()
