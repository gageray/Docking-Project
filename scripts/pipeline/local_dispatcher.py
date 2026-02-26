import json
import os
import itertools
from pathlib import Path



import argparse

def poll_google_drive(creds_path: str) -> set[str]:
    """
    Connect to Google Drive and return a set of completed output filenames.
    """
    import sys
    sys.path.append(str(Path("scripts/kaggle").absolute()))
    import drive_auth, drive_io # type: ignore
    
    # Needs a local credentials file
    if not os.path.exists(creds_path):
        print(f"[-] WARNING: Local Google Drive credentials not found at {creds_path}")
        print("[-] Skipping Drive polling. Assuming 0 completed jobs.")
        return set()
        
    try:
        service, _ = drive_auth.setup_drive(creds_path, verify=False)
        files = drive_io.list_drive_folder(service, drive_auth.DEFAULT_FOLDERS["outputs"], verbose=False)
        # GNINA outputs are named like 6X3X_apo_mol_001_out.sdf
        return set(f['name'] for f in files)
    except Exception as e:
        print(f"[-] Error polling Google Drive: {e}")
        return set()

def load_box_metadata(box_config_path: Path) -> dict:
    """Load the manual box boundaries."""
    if not box_config_path.exists():
        return {"__default__": {"center_x": 0.0, "center_y": 0.0, "center_z": 0.0, "size_x": 25.0, "size_y": 25.0, "size_z": 25.0}}
    with open(box_config_path, 'r') as f:
        data = json.load(f)
        return data if isinstance(data, dict) else {}

def discover_local_files(receptor_dir: Path, ligand_dir: Path):
    """Discover prepped receptors and ligands ready for docking."""
    receptors = [f.name for f in receptor_dir.glob("*_apo.pdbqt")] if receptor_dir.exists() else []
    ligands = [f.name for f in ligand_dir.glob("*.sdf")] if ligand_dir.exists() else []
    return receptors, ligands

def build_job_queue(box_config_path: Path, receptor_dir: Path, ligand_dir: Path, creds_path: str):
    """
    Reads local configurations and generates a fully resolved array of GNINA parameters.
    """
    box_cfg = load_box_metadata(box_config_path)
    receptors, ligands = discover_local_files(receptor_dir, ligand_dir)
    completed = poll_google_drive(creds_path)
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

def upload_queue_to_drive(work_queue_path: Path, creds_path: str):
    """Upload work_queue.json to Google Drive."""
    import sys
    sys.path.append(str(Path("scripts/kaggle").absolute()))
    import drive_auth, drive_io # type: ignore
    
    if not os.path.exists(creds_path):
        print(f"[-] FATAL: Local Google Drive credentials not found at {creds_path}")
        print("[-] Cannot upload work queue.")
        sys.exit(1)
        
    try:
        service, _ = drive_auth.setup_drive(creds_path, verify=False)
        print(f"[*] Uploading {work_queue_path} to Google Drive...")
        drive_io.upload_to_drive(service, str(work_queue_path), drive_auth.DEFAULT_FOLDERS["root"], overwrite=True)
    except Exception as e:
        print(f"[-] Error uploading to Google Drive: {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="GNINA Kaggle Worker Dispatcher")
    parser.add_argument("--config", type=Path, default=Path("scripts/pipeline/local_dispatcher.json"), help="Path to main pipeline JSON configuration")
    parser.add_argument("--box-config", type=Path, help="Path to box configuration JSON")
    parser.add_argument("--receptor-dir", type=Path, help="Directory containing prepped receptors")
    parser.add_argument("--ligand-dir", type=Path, help="Directory containing ligand SDFs")
    parser.add_argument("--work-queue", type=Path, help="Output path for the generated work queue")
    parser.add_argument("--creds", type=str, help="Path to Google Drive API credentials")
    parser.add_argument("--kaggle-dir", type=Path, help="Directory containing Kaggle Kernel scripts and metadata")
    
    args = parser.parse_args()
    
    config_data = {}
    if args.config and args.config.exists():
        with open(args.config, 'r') as f:
            config_data = json.load(f)
            
    # Priority: CLI Argument > JSON Config > Hardcoded Default
    box_config = args.box_config or Path(config_data.get("box_config", "scripts/scheduler/box_config.json"))
    receptor_dir = args.receptor_dir or Path(config_data.get("receptor_dir", "data/receptors/prepped"))
    ligand_dir = args.ligand_dir or Path(config_data.get("ligand_dir", "data/ligands"))
    work_queue = args.work_queue or Path(config_data.get("work_queue", "work_queue.json"))
    creds = args.creds or config_data.get("creds", "credentials/key.json")
    kaggle_dir = args.kaggle_dir or Path(config_data.get("kaggle_dir", "scripts/kaggle"))
    
    print("Checking for incomplete jobs...")
    jobs = build_job_queue(box_config, receptor_dir, ligand_dir, creds)
    
    if not jobs:
        print("All docking jobs completed! Nothing to push to Kaggle.")
        return
        
    print(f"Submitting {len(jobs)} completely resolved jobs to Kaggle cluster...")
    
    with open(work_queue, 'w') as f:
        json.dump(jobs, f, indent=2)
        
    print(f"Wrote {work_queue}.")
    
    upload_queue_to_drive(work_queue, creds)
    
    print("[*] Pushing kernel to Kaggle...")
    # NOTE: user must have kaggle installed locally and configured via ~/.kaggle/kaggle.json
    import subprocess
    try:
        subprocess.run(["kaggle", "kernels", "push", "-p", str(kaggle_dir)], check=True)
        print("[+] Kernel successfully pushed to Kaggle!")
    except Exception as e:
        print(f"[-] Failed to push kernel: {e}")
    
if __name__ == "__main__":
    main()
