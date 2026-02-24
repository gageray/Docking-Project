import os
import json
import pymol
from pymol import cmd
import math

def process_receptor(cif_path, meta_path, out_dir):
    if not os.path.exists(cif_path) or not os.path.exists(meta_path):
        print(f"Skipping {cif_path} (missing file or metadata).")
        return

    print(f"Processing {cif_path}...")
    with open(meta_path, 'r') as f:
        meta = json.load(f)
        
    cmd.load(cif_path, "complex")
    
    target_lig_resn = meta.get("target_ligand_resn")
    target_lig_chain_out = meta.get("target_ligand_chain") # e.g. 'Z'
    
    # Identify the receptor chains to keep
    chain_strs = []
    for c in meta.get("chains", {}).keys():
        chain_strs.append(f"chain {c}")
    receptor_sel = "(" + " or ".join(chain_strs) + ")"
    
    if target_lig_resn and target_lig_chain_out:
        # 1. Identify BZD chains to compute center
        bzd_chains = [c for c, desc in meta.get("chains", {}).items() if "_bzd" in desc]
        if not bzd_chains:
            print("Warning: No _bzd chains found in metadata, cannot compute center.")
            cmd.delete("all")
            return

        bzd_sel = " or ".join([f"chain {c}" for c in bzd_chains])
        bzd_com = cmd.centerofmass(bzd_sel)
        
        # 2. Find closest ligand to the compute center
        model = cmd.get_model(f"resn {target_lig_resn}")
        ligands = set((a.chain, a.resi) for a in model.atom)
        
        closest_lig_sel = None
        min_dist = float('inf')
        
        for ch, resi in ligands:
            lig_sel = f"(chain {ch} and resi {resi} and resn {target_lig_resn})"
            lig_com = cmd.centerofmass(lig_sel)
            dist = math.dist(bzd_com, lig_com)
            if dist < min_dist:
                min_dist = dist
                closest_lig_sel = lig_sel
                
        if not closest_lig_sel:
            print(f"Warning: Could not find any {target_lig_resn} in structure.")
            cmd.delete("all")
            return
            
        print(f"Identified closest BZD ligand: {closest_lig_sel} at distance {min_dist:.2f}")
        
        # 3. Alter its chain to Z (or whatever target_ligand_chain is)
        cmd.alter(closest_lig_sel, f"chain='{target_lig_chain_out}'")
        cmd.sort()
        
        # 4. Prune everything but the target receptor chains and the target ligand (now chain Z)
        final_lig_sel = f"(resn {target_lig_resn} and chain {target_lig_chain_out})"
        cmd.remove(f"not ({receptor_sel} or {final_lig_sel})")
        # Ensure no other residue with this resn exists (duplicates)
        cmd.remove(f"resn {target_lig_resn} and not {final_lig_sel}")
        
        out_name = os.path.basename(cif_path).split('.')[0].replace(" (legacy diazepam)", "")
        out_pdb = os.path.join(out_dir, f"{out_name}_ligand.pdb")
    else:
        cmd.remove(f"not {receptor_sel}")
        out_name = os.path.basename(cif_path).split('.')[0]
        out_pdb = os.path.join(out_dir, f"{out_name}_apo.pdb")
        
    cmd.save(out_pdb, "complex")
    print(f"Saved {out_pdb}")
    cmd.delete("all")

def main():
    pymol.pymol_argv = ['pymol', '-c', '-q']
    pymol.finish_launching()
    
    base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    receptors_dir = os.path.join(base_dir, "data", "receptors")
    
    receptors = [
        ("6X3U.cif", "6x3u_metadata.json"),
        ("6X3X (legacy diazepam).cif", "6x3x_metadata.json"),
        ("9CRS.cif", "9crs_metadata.json"),
        ("9CSB.cif", "9csb_metadata.json")
    ]
    
    for cif_name, meta_name in receptors:
        process_receptor(
            os.path.join(receptors_dir, cif_name),
            os.path.join(receptors_dir, meta_name),
            receptors_dir
        )
        
    cmd.quit()

if __name__ == "__main__":
    main()
