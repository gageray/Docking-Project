"""
Post-Docking Pose Analysis
Run inside PyMOL: `run scripts/vis_results.py`
"""
from pymol import cmd
import glob
import os
import json

def load_local_config(config_path="config.json"):
    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            return json.load(f)
    return {}

def analyze_docking_results(target_pdb=None, docking_results_path="data/docking_results/*.sdf"):
    cmd.reinitialize()
    config = load_local_config()
    
    # 1. Provide an option to hardcode or glob the receptor and results
    # For now, we assume the user might have to specify or we load a default like 9CSB
    if target_pdb is None:
        target_pdb = os.path.join("data", "receptors", "9CSB_prepped.pdb")
        
    if not os.path.exists(target_pdb):
        print(f"Receptor file {target_pdb} not found. Ensure you have run prep.")
        return
        
    cmd.load(target_pdb, "receptor")
    print(f"Loaded receptor: {target_pdb}")
    
    # 2. Setup receptor representation
    cmd.hide("everything", "receptor")
    cmd.show("cartoon", "receptor")
    cmd.color("gray70", "receptor")
    
    # Show surface for clash detection
    cmd.show("surface", "receptor")
    cmd.set("transparency", 0.3)
    
    # Highlight Alpha1-His102 (Requires proper chain/resi assignment based on your PDB)
    # Often His102 is exactly resi 102, but check your specific structural numbering.
    his_resi = config.get("visualization", {}).get("key_residue", "102")
    cmd.select("key_his", f"receptor and resn HIS and resi {his_resi}")
    if cmd.count_atoms("key_his") > 0:
        cmd.show("sticks", "key_his")
        cmd.color("magenta", "key_his")
        print("Highlighted His102")
        
    # 3. Load top docked poses
    files = glob.glob(docking_results_path)
    if not files:
        print(f"No docking results found at {docking_results_path}")
        print("Visualization script is ready for when GNINA outputs exist.")
        return
        
    print(f"Loading {len(files)} docked output files...")
    
    for i, f in enumerate(files):
        obj_name = f"pose_{i+1}"
        cmd.load(f, obj_name)
        cmd.show("sticks", obj_name)
        cmd.color("green", obj_name)
        
        # 4. Measure generic distances (Example: Carbonyl Oxygen to Receptor)
        # Assuming the ester oxygen is identifiable. If GNINA outputs retain atom order,
        # or we just find the closest O to the receptor.
        # This is a broad selection for demonstration: finding H-bonds between ligand and receptor
        cmd.distance(f"hbond_{i+1}", f"{obj_name} and name O*", "receptor and (name N* or name O*)", cutoff=3.2, mode=2)
        
        # 5. Highlight severe clashes (Ligand atoms within 1.5A of receptor atoms that are not H-bonding)
        clash_name = f"clash_{i+1}"
        cmd.select(clash_name, f"{obj_name} within 1.5 of receptor")
        if cmd.count_atoms(clash_name) > 0:
            cmd.color("red", clash_name)
            cmd.show("spheres", clash_name)
            cmd.set("sphere_scale", 0.5, clash_name)
            print(f"Warning: {cmd.count_atoms(clash_name)} clashing atoms detected in {obj_name}")
            
    cmd.zoom("pose_*", buffer=10.0)
    print("Success Criteria Check: Ferulic acid headgroup should sit deep in the pocket, and C22 tail project into solvent.")

if __name__ in ["__main__", "pymol"]:
    analyze_docking_results()
