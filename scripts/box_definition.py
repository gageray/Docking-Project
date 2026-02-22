import json
import os
import sys

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, Superimposer
from Bio import Align
from Bio.PDB.Polypeptide import three_to_one

from utils import load_config, setup_logger

logger = setup_logger(__name__)

def get_ligand_centroid(structure, resname):
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname == resname:
                    for atom in residue:
                        coords.append(atom.get_coord())
    if not coords:
        return None
    return np.mean(coords, axis=0)

def extract_sequence_and_cas(structure, chain_ids):
    seq_cas_pairs = []
    for model in structure:
        for chain in model:
            if chain.id in chain_ids:
                for residue in chain:
                    if "CA" in residue:
                        try:
                            # Try to get 1-letter code
                            res_code = three_to_one(residue.resname)
                            seq_cas_pairs.append((res_code, residue["CA"]))
                        except KeyError:
                            # Non-standard amino acid without 3-to-1 mapping
                            continue
    return seq_cas_pairs

def main():
    try:
        config = load_config()
    except Exception as e:
        logger.info(f"Failed to load config.json: {e}")
        return
        
    box_cfg = config.get("box_definition", {})
    if not box_cfg:
        logger.info("No box_definition config found.")
        return
        
    ref_cfg = box_cfg.get("reference_box", {})
    ref_cif = ref_cfg.get("raw_cif")
    ligand_resname = ref_cfg.get("ligand_resname")
    ref_pdb = ref_cfg.get("prepped_pdb")
    ref_chains = ref_cfg.get("pocket_chains", [])
    
    out_file = box_cfg.get("output_file", "box_params.json")
    
    # 1. Reference Centroid
    cif_parser = MMCIFParser(QUIET=True)
    try:
        ref_raw_struct = cif_parser.get_structure("ref_raw", ref_cif)
        ref_centroid = get_ligand_centroid(ref_raw_struct, ligand_resname)
        if ref_centroid is None:
            logger.info(f"Error: Could not find ligand {ligand_resname} in {ref_cif}")
            return
        logger.info(f"Reference centroid ({ligand_resname}): {ref_centroid}")
    except Exception as e:
        logger.info(f"Failed to calculate reference centroid: {e}")
        return
        
    # Load prepared reference structure for sequence/CA extraction
    pdb_parser = PDBParser(QUIET=True)
    try:
        ref_prepped_struct = pdb_parser.get_structure("ref_prepped", ref_pdb)
        ref_seq_cas = extract_sequence_and_cas(ref_prepped_struct, ref_chains)
        ref_seq = "".join([x[0] for x in ref_seq_cas])
        ref_cas = [x[1] for x in ref_seq_cas]
        if not ref_seq:
            logger.info("Error: Could not extract sequence from reference PDB.")
            return
    except Exception as e:
        logger.info(f"Failed to parse or extract seq from reference PDB {ref_pdb}: {e}")
        return
        
    aligner = Align.PairwiseAligner()
    box_params = {
        "reference": {
            "name": "reference",
            "centroid": ref_centroid.tolist()
        },
        "targets": []
    }
    
    targets = box_cfg.get("targets", [])
    
    for target in targets:
        t_name = target.get("name")
        t_pdb = target.get("prepped_pdb")
        t_chains = target.get("pocket_chains", [])
        
        logger.info(f"\nAligning target {t_name} from {t_pdb}")
        
        if not os.path.exists(t_pdb):
            logger.info(f"Target PDB {t_pdb} does not exist. Skipping.")
            continue
            
        try:
            t_struct = pdb_parser.get_structure(t_name, t_pdb)
            t_seq_cas = extract_sequence_and_cas(t_struct, t_chains)
            t_seq = "".join([x[0] for x in t_seq_cas])
            t_cas = [x[1] for x in t_seq_cas]
            
            if not t_seq:
                logger.info(f"No sequence found for target {t_name}")
                continue
                
            alignment = aligner.align(ref_seq, t_seq)[0]
            
            ref_aligned_cas = []
            target_aligned_cas = []
            
            # Map indices from alignment blocks
            for ref_block, t_block in zip(alignment.aligned[0], alignment.aligned[1]):
                for i, j in zip(range(ref_block[0], ref_block[1]), range(t_block[0], t_block[1])):
                    ref_aligned_cas.append(ref_cas[i])
                    target_aligned_cas.append(t_cas[j])
                    
            if len(ref_aligned_cas) < 3:
                logger.info(f"Error: Not enough aligned atoms ({len(ref_aligned_cas)}) for target {t_name}")
                continue
                
            # Superimposer target onto ref to get translation/rotation matrix
            # Spec states: "Apply the inverse translation/rotation matrix to the reference centroid to map it into the target space"
            # Which implies ref is fixed, target is moving, and then applying the inverse.
            # Here, we swap the argument order: TARGET is fixed, REF is moving, and we apply rotran forward.
            # This is mathematically equivalent and avoids calculating the inverse matrix.
            super_imposer = Superimposer()
            super_imposer.set_atoms(target_aligned_cas, ref_aligned_cas)
            
            rmsd = super_imposer.rms
            if rmsd > 3.0:
                logger.info(f"WARNING: RMSD for {t_name} is high: {rmsd:.2f} A")
            else:
                logger.info(f"RMSD: {rmsd:.2f} A")
                
            rot, tran = super_imposer.rotran
            # Mapping ref centroid to target:
            # rot is a Numpy matrix. The mapping is: new_coord = np.dot(coord, rot) + tran
            t_centroid = np.dot(ref_centroid, rot) + tran
            t_centroid = np.array(t_centroid, dtype=float)
            
            logger.info(f"Mapped target centroid: {t_centroid}")
            
            box_params["targets"].append({
                "name": t_name,
                "rmsd": float(rmsd),
                "centroid": t_centroid.tolist()
            })
            
        except Exception as e:
            logger.info(f"Failed alignment for target {t_name}: {e}")
            continue
            
    # Export params
    try:
        with open(out_file, "w") as f:
            json.dump(box_params, f, indent=4)
        logger.info(f"\nSuccessfully wrote parameters to {out_file}")
    except Exception as e:
        logger.info(f"Failed to write output {out_file}: {e}")

if __name__ == "__main__":
    main()
