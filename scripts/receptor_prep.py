import json
import os
import sys
import tempfile

from Bio.PDB import MMCIFParser, PDBIO, Select

try:
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False
    logger.info("Warning: openmm or pdbfixer is not installed. receptor_prep requires these to run completely.")

class ChainSelector(Select):
    def __init__(self, keep_chains):
        self.keep_chains = keep_chains

    def accept_chain(self, chain):
        if chain.get_id() in self.keep_chains:
            return 1
        return 0

from utils import load_config, setup_logger

logger = setup_logger(__name__)

def process_receptor(receptor_item):
    in_cif = receptor_item.get("input")
    out_pdb = receptor_item.get("output")
    keep_chains = receptor_item.get("keep_chains", [])
    
    if not os.path.exists(in_cif):
        logger.info(f"Error: Target input {in_cif} does not exist.")
        return

    logger.info(f"Processing {in_cif} -> {out_pdb}")
    
    # 1. Load CIF via Biopython MMCIFParser
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure("struct", in_cif)
    except Exception as e:
        logger.info(f"Failed to parse CIF {in_cif}: {e}")
        return

    # 2. Delete all chains not explicitly listed and map to intermediate PDB
    io = PDBIO()
    io.set_structure(structure)
    
    selector = ChainSelector(keep_chains)
    
    # 3. Save intermediate PDB to temp file
    temp_pdb = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    temp_pdb.close()
    
    try:
        io.save(temp_pdb.name, select=selector)
        
        # 4. Run PDBFixer
        if not HAS_OPENMM:
            raise ImportError("openmm and pdbfixer are required to process receptors.")
        fixer = PDBFixer(filename=temp_pdb.name)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.4)
        fixer.removeHeterogens(keepWater=False)
        
        # 5. Write final PDB to output path
        os.makedirs(os.path.dirname(out_pdb), exist_ok=True)
        with open(out_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
            
        # Log output atom and chain counts
        num_atoms = sum(1 for _ in fixer.topology.atoms())
        num_chains = sum(1 for _ in fixer.topology.chains())
        logger.info(f"[{out_pdb}] Success. Atoms: {num_atoms}, Chains: {num_chains}")
        
    except Exception as e:
        logger.info(f"Failed to fix/export {in_cif}: {e}")
    finally:
        os.unlink(temp_pdb.name)

def main():
    try:
        config = load_config()
    except Exception as e:
        logger.info(f"Failed to load config.json: {e}")
        return
        
    receptor_list = config.get("receptor_prep", [])
    
    for item in receptor_list:
        try:
            process_receptor(item)
        except Exception as e:
            logger.info(f"Unexpected error processing item {item}: {e}")
            continue

if __name__ == "__main__":
    main()
