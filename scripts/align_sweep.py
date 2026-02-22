import glob
import os
import csv
import logging
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Geometry import Point3D
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

def get_centroid(conf, atom_indices):
    coords = [conf.GetAtomPosition(idx) for idx in atom_indices]
    coords = np.array([[pt.x, pt.y, pt.z] for pt in coords])
    centroid = np.mean(coords, axis=0)
    return Point3D(*centroid)

def main():
    data_dir = os.path.join(os.path.dirname(__file__), "..", "data", "ligands")
    out_sdf = os.path.join(data_dir, "df_aligned_sweep.sdf")
    out_csv = os.path.join(data_dir, "df_sweep_reach.csv")
    
    # 1. Load Conformers
    sdf_files = glob.glob(os.path.join(data_dir, "docosanyl_ferulate_conf*.sdf"))
    
    # Sort files numerically by number at the end
    def extract_conf_num(filepath):
        basename = os.path.basename(filepath)
        num_str = basename.replace("docosanyl_ferulate_conf", "").replace(".sdf", "")
        return int(num_str) if num_str.isdigit() else 999
        
    sdf_files.sort(key=extract_conf_num)
    
    if not sdf_files:
        logging.error("No conformer files found!")
        return
        
    mols = []
    for f in sdf_files:
        suppl = Chem.SDMolSupplier(f)
        for m in suppl:
            if m is not None:
                m.SetProp("_Name", os.path.basename(f))
                mols.append(m)
                
    if not mols:
        logging.error("No valid molecules loaded.")
        return
        
    # 2. Map the Pharmacophore (Rigid Phenol Ring Anchor)
    # The C=C-C=O linkage was varied in torsion during the grid search.
    # Therefore, we must align ONLY on the absolutely rigid part: the phenol ring and the first alkene carbon.
    # SMARTS for the Guaiacol ring + 1st alkene carbon: c1cc(O)c(OC)cc1C
    fa_core_smarts = 'c1cc(O)c(OC)cc1C'
    core_q = Chem.MolFromSmarts(fa_core_smarts)
        
    match_indices = mols[0].GetSubstructMatch(core_q)
    if not match_indices:
        logging.error(f"Could not match rigid phenol anchor in first molecule!")
        return

    logging.info(f"Matched FA Core to {len(match_indices)} atoms.")
    atom_map = [(idx, idx) for idx in match_indices]
    
    # Let's also find the C22 terminal carbon. 
    # It is a primary aliphatic CH3 attached to a CH2. 'CC' where both are sp3.
    # Better: '[CX4H3][CX4H2]'
    terminal_q = Chem.MolFromSmarts('[CX4H3][CX4H2]')
    term_matches = mols[0].GetSubstructMatches(terminal_q)
    c22_idx = None
    if term_matches:
        # We might match the methoxy group if we aren't careful, but methoxy is [CX4H3] attached to Oxygen, not [CX4H2].
        # So [CX4H3][CX4H2] will only match the docosanyl tail.
        # Check which of the two matched atoms is the terminal one (the CX4H3).
        for match in term_matches:
            atom = mols[0].GetAtomWithIdx(match[0])
            if atom.GetDegree() == 1:
                c22_idx = match[0]
                break
    
    if c22_idx is None:
        logging.warning("Could not find terminal C22 carbon definitively. Distance calculations might fail.")

    ref_mol = mols[0]
    
    writer = Chem.SDWriter(out_sdf)
    csv_f = open(out_csv, 'w', newline='')
    csv_writer = csv.writer(csv_f)
    csv_writer.writerow(['Conformer', 'C22_Distance_Angstroms'])
    
    # 3 & 4: Coordinate Alignment & Verification
    for i, mol in enumerate(mols):
        if i == 0:
            rmsd = 0.0
        else:
            rmsd = rdMolAlign.AlignMol(mol, ref_mol, atomMap=atom_map)
            
        # Verify - increased tolerance because MMFF optimization allows for minor coordinate drift (<0.05 A)
        if rmsd > 0.05 and i > 0:
            logging.warning(f"Mol {mol.GetProp('_Name')}: RMSD of FA rigid anchor is {rmsd:.4f} > 0.05 Å")
            
        mname = mol.GetProp('_Name') if mol.HasProp('_Name') else f"conf_{i}"
        
        # Distance calculation
        dist = "N/A"
        if c22_idx is not None:
            # centroid of FA core in the newly aligned molecule
            conf = mol.GetConformer()
            centroid = get_centroid(conf, match_indices)
            c22_pos = conf.GetAtomPosition(c22_idx)
            dist = centroid.Distance(c22_pos)
            csv_writer.writerow([mname, f"{dist:.4f}"])
        else:
            csv_writer.writerow([mname, "N/A"])
            
        writer.write(mol)
        
    writer.close()
    csv_f.close()
    
    logging.info(f"Done! Wrote {len(mols)} aligned conformers to {out_sdf} and metrics to {out_csv}")
    
    print("\n--- PyMOL Rendering Script Snippet ---")
    print("# To be executed in PyMOL after loading df_aligned_sweep.sdf")
    print("cmd.set('all_states', 1)")
    print("cmd.hide('everything')")
    print("cmd.show('sticks', 'df_aligned_sweep')")
    print("cmd.color('atomic', 'df_aligned_sweep')")
    print(f"cmd.select('tail', 'not id {'+'.join(str(i+1) for i in match_indices)}') # PyMOL uses 1-based indexing for limits")
    print("cmd.color('gray50', 'tail')")
    print("cmd.set('stick_transparency', 0.7, 'tail')")

if __name__ == '__main__':
    main()
