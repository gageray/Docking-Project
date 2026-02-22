import json
import os
import sys
import itertools
import math
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    import dimorphite_dl
except ImportError:
    dimorphite_dl = None

from utils import load_config, setup_logger

logger = setup_logger(__name__)

def protonate(smiles):
    if dimorphite_dl is None:
        logger.info("DimorphiteDL not installed, skipping protonation.")
        return smiles
    
    try:
        # dimorphite_dl.protonate_smiles returns a list/iterator of smiles
        protonated = list(dimorphite_dl.protonate_smiles(
            smiles, 
            ph_min=7.4, 
            ph_max=7.4, 
            precision=0.0
        ))
        
        if protonated:
            logger.info(f"Protonated {smiles} -> {protonated[0]}")
            return protonated[0]
    except Exception as e:
        logger.info(f"Failed to protonate {smiles}: {e}")
    return smiles

def set_dihedral(conf, atom_indices, angle_deg):
    i, j, k, l = atom_indices
    Chem.rdMolTransforms.SetDihedralDeg(conf, i, j, k, l, angle_deg)

def process_single(name, smiles, outdir):
    prot_smiles = protonate(smiles)
    mol = Chem.MolFromSmiles(prot_smiles)
    if mol is None:
        logger.info(f"Invalid SMILES for {name}: {prot_smiles}")
        return
    
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    AllChem.EmbedMolecule(mol, params)
    
    if AllChem.MMFFOptimizeMolecule(mol, maxIters=2000) != 0:
        logger.info(f"Warning: MMFF optimization failed to converge for {name}")
    
    outpath = os.path.join(outdir, f"{name}.sdf")
    writer = Chem.SDWriter(outpath)
    writer.write(mol)
    writer.close()
    logger.info(f"[{name}] Exported to {outpath}")

    # ester_match = (c_beta, c_alpha, c_carbonyl, o_carbonyl, o_ester, c_alkyl1)
    # in the SMARTS [CX3H1]=[CX3H1]C(=O)[OX2][CH2], indices are:
    # 0: c_beta, 1: c_alpha, 2: c_carbonyl, 3: o_carbonyl, 4: o_ester, 5: c_alkyl1
    # But wait, our previous code assumed: (c_carbonyl, o_carbonyl, o_ester, c_alkyl1)
    # The actual passed match tuple is 6 elements: (c_beta, c_alpha, c_carbonyl, o_carbonyl, o_ester, c_alkyl1)
    # Let's just use the indices directly from whatever is passed in, assuming the last 4 are the ester core.
    
def identify_tail_atoms(mol, ester_match):
    # Let's extract based on the signature length (6 from script, 6 from test)
    c_carbonyl_idx = ester_match[-4]
    o_ester_idx = ester_match[-2]
    c_alkyl1_idx = ester_match[-1]
    
    tail_atoms = []
    
    # We pass visited set to prevent loops
    def dfs(curr_idx, visited):
        visited.add(curr_idx)
        tail_atoms.append(curr_idx)
        atom = mol.GetAtomWithIdx(curr_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited and neighbor.GetAtomicNum() == 6:  # only follow carbons
                dfs(n_idx, set(visited))
    
    # Start DFS at the first carbon of the tail, block traversal back to the ester oxygen and carbonyl
    dfs(c_alkyl1_idx, {o_ester_idx, c_carbonyl_idx})
    return tail_atoms

def process_grid_search(name, smiles, outdir):
    prot_smiles = protonate(smiles)
    mol = Chem.MolFromSmiles(prot_smiles)
    if mol is None: return
    mol = Chem.AddHs(mol)
    
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    
    # Identify linkage using SMARTS C(=O)O[CH2] -> more robust: [CX3](=[OX1])[OX2][CH2X4]
    # But let's find the specific dihedrals needed for sweep (from spec):
    # 1. C=C-C=O (sp2-sp2)
    # 2. O-CH2 (sp3-sp3)
    # 3. CH2-CH2 (sp3-sp3)
    
    # Spec: "Match ester linkage C(=O)O[CH2]"
    # Use robust SMARTS to avoid matching aromatic rings
    patt = Chem.MolFromSmarts('[CX3H1]=[CX3H1]C(=O)[OX2][CH2]')
    match = mol.GetSubstructMatch(patt)
    
    if not match:
        logger.info(f"[{name}] Could not find matching linkage for grid search. Falling back to single.")
        process_single(name, smiles, outdir)
        return
        
    c_beta, c_alpha, c_carbonyl, o_carbonyl, o_ester, c_alkyl1 = match[:6]
    
    # Identify tail to keep straight
    tail_carbons = identify_tail_atoms(mol, (c_carbonyl, o_carbonyl, o_ester, c_alkyl1))
    
    # Force tail to all trans (180 deg)
    conf = mol.GetConformer()
    for i in range(len(tail_carbons) - 3):
        idx1, idx2, idx3, idx4 = tail_carbons[i:i+4]
        set_dihedral(conf, (idx1, idx2, idx3, idx4), 180.0)
    
    # Define sweep pivots
    # 1. C=C-C=O  (c_beta - c_alpha - c_carbonyl - o_ester)
    # 2. C-O-C-C  (c_carbonyl - o_ester - c_alkyl1 - c_alkyl2)  (O-CH2)
    # 3. O-C-C-C  (o_ester - c_alkyl1 - c_alkyl2 - c_alkyl3)    (CH2-CH2)
    
    if len(tail_carbons) < 3:
        logger.info(f"[{name}] Tail too short for sweep. Falling back.")
        process_single(name, smiles, outdir)
        return
        
    c_alkyl2 = tail_carbons[1]
    c_alkyl3 = tail_carbons[2] if len(tail_carbons) > 2 else None
    
    angles_1 = [0, 180]
    angles_2 = [0, 60, 120, 180, 240, 300]
    angles_3 = [0, 60, 120, 180, 240, 300]
    
    combinations = list(itertools.product(angles_1, angles_2, angles_3))
    
    logger.info(f"[{name}] Generating {len(combinations)} conformers via grid search...")
    
    results = []
    
    for idx, (a1, a2, a3) in enumerate(combinations):
        # Create a new conformer for each combination
        new_mol = Chem.Mol(mol)
        new_conf = new_mol.GetConformer()
        
        # Apply the dihedrals
        set_dihedral(new_conf, (c_beta, c_alpha, c_carbonyl, o_ester), a1)
        set_dihedral(new_conf, (c_carbonyl, o_ester, c_alkyl1, c_alkyl2), a2)
        if c_alkyl3 is not None:
            set_dihedral(new_conf, (o_ester, c_alkyl1, c_alkyl2, c_alkyl3), a3)
            
        # Optimize with fixed head/tail
        props = AllChem.MMFFGetMoleculeProperties(new_mol)
        if props is None: continue
        ff = AllChem.MMFFGetMoleculeForceField(new_mol, props, confId=0)
        
        # Add constraints for tail carbons to stay rigidly straight
        for c in tail_carbons:
            # maxDispl, forceConstant: keeping force constant high but displacement tight
            ff.MMFFAddPositionConstraint(c, 0.0, 1.0e4)
            
        ff.Initialize()
        converged = ff.Minimize(maxIts=2000)
        
        if converged == 0:
            energy = ff.CalcEnergy()
            results.append((energy, new_mol, idx))
            
    if not results:
        logger.info(f"[{name}] All conformers failed to converge.")
        return
        
    results.sort(key=lambda x: x[0])
    global_min = results[0][0]
    
    valid_results = [res for res in results if res[0] <= global_min + 20.0]
    
    for i, (energy, conf_mol, orig_idx) in enumerate(valid_results):
        outpath = os.path.join(outdir, f"{name}_conf{i+1}.sdf")
        writer = Chem.SDWriter(outpath)
        writer.write(conf_mol)
        writer.close()
        
    logger.info(f"[{name}] Exported {len(valid_results)} conformers (dropped {len(combinations) - len(valid_results)}).")

def main():
    try:
        config = load_config()
    except Exception as e:
        logger.info(f"Failed to load config.json: {e}")
        return
        
    ligand_cfg = config.get("ligand_prep", {})
    outdir = ligand_cfg.get("output_dir", "ligands")
    os.makedirs(outdir, exist_ok=True)
    
    ligands = ligand_cfg.get("ligands", {})
    
    for name, data in ligands.items():
        smiles = data.get("smiles")
        mode = data.get("mode", "single")
        
        try:
            if mode == "single":
                process_single(name, smiles, outdir)
            elif mode == "grid_search":
                process_grid_search(name, smiles, outdir)
            else:
                logger.info(f"Unknown mode '{mode}' for {name}. Using single.")
                process_single(name, smiles, outdir)
        except Exception as e:
            logger.info(f"Error processing {name}: {e}")
            continue

if __name__ == "__main__":
    main()
