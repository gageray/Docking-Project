import json
import os
import itertools
from pathlib import Path

# --- Configuration ---
BOX_CONFIG_PATH = Path("scripts/scheduler/box_config.json")
RECEPTOR_DIR = Path("data/receptors/prepped")
LIGAND_DIR = Path("data/ligands")
WORK_QUEUE_PATH = Path("work_queue.json")
KERNEL_METADATA_PATH = Path("kernel-metadata.json")

def poll_google_drive():
    """
    Connect to Google Drive and return a set of completed output filenames.
    Currently a mock holding space until we copy the drive_auth logic.
    """
    print("Polling Google Drive for completed jobs...")
    return set()

def load_box_metadata() -> dict:
    """Load the manual box boundaries."""
    if not BOX_CONFIG_PATH.exists():
        return {"__default__": {"center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": 25.0, "size_y": 25.0, "size_z": 25.0}}
    with open(BOX_CONFIG_PATH, 'r') as f:
        data = json.load(f)
        return data if isinstance(data, dict) else {}

def discover_local_files():
    """Discover prepped receptors and ligands ready for docking."""
    receptors = [f.name for f in RECEPTOR_DIR.glob("*_apo.pdbqt")] if RECEPTOR_DIR.exists() else ["6X3X_apo.pdbqt", "4LDE_apo.pdbqt"]
    ligands = [f.name for f in LIGAND_DIR.glob("*.sdf")] if LIGAND_DIR.exists() else ["mol_001.sdf", "mol_002.sdf"]
    return receptors, ligands

def build_job_queue():
    """
    Reads local configurations and generates a fully resolved array of GNINA parameters.
    """
    box_cfg = load_box_metadata()
    receptors, ligands = discover_local_files()
    completed = poll_google_drive()
    queue = []
    
    # Ensure __default__ exists for the Pyre linter
    default_box = box_cfg.get("__default__", {})
    if not isinstance(default_box, dict):
        default_box = {}
        
    for rec, lig in itertools.product(receptors, ligands):
        out_name = f"{rec.replace('.pdbqt', '')}_{lig.replace('.sdf', '')}_out.sdf"
        
        if out_name in completed:
            continue
            
        rec_key = rec.replace('_apo.pdbqt', '_prepped')
        
        box = box_cfg.get(rec_key, default_box)
        if not isinstance(box, dict):
            box = default_box
            
        queue.append({
            "receptor": rec,
            "ligand": lig,
            "cx": float(box.get("center_x", 0.0)),
            "cy": float(box.get("center_y", 0.0)),
            "cz": float(box.get("center_z", 0.0)),
            "sx": float(box.get("size_x", 25.0)),
            "sy": float(box.get("size_y", 25.0)),
            "sz": float(box.get("size_z", 25.0)),
        })
        
        # Respect Kaggle batch limits (e.g. 500 per push)
        if len(queue) >= 500:
            break
            
    return queue

def generate_kaggle_metadata():
    """Generate the kernel-metadata.json required for the API push."""
    metadata = {
      "id": "your-username/gnina-docking-job",
      "title": "GNINA Docking Job",
      "code_file": "gnina_worker.py",
      "language": "python",
      "kernel_type": "script",
      "is_private": "true",
      "enable_gpu": "true",
      "enable_internet": "true",
      "dataset_sources": [],
      "competition_sources": [],
      "kernel_sources": [],
      "model_sources": []
    }
    with open(KERNEL_METADATA_PATH, 'w') as f:
        json.dump(metadata, f, indent=2)

def main():
    print("Checking for incomplete jobs...")
    jobs = build_job_queue()
    
    if not jobs:
        print("All docking jobs completed! Nothing to push to Kaggle.")
        return
        
    print(f"Submitting {len(jobs)} completely resolved jobs to Kaggle cluster...")
    
    with open(WORK_QUEUE_PATH, 'w') as f:
        json.dump(jobs, f, indent=2)
        
    generate_kaggle_metadata()
    print(f"Wrote {WORK_QUEUE_PATH} and {KERNEL_METADATA_PATH}.")
    print("Next step: cd into the push directory and run 'kaggle kernels push'")
    
if __name__ == "__main__":
    main()
