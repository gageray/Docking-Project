#!/usr/bin/env python3
"""Upload files to Drive push folder for Kaggle docking jobs.

Usage:
    python build_kaggle_push.py --receptor 6X3U_apo.pdbqt --ligand flumazenil.sdf --box 0,0,0,20,20,20
    python build_kaggle_push.py --config job_config.json
"""

import json
import shutil
import argparse
import sys
from pathlib import Path
from googleapiclient.http import MediaFileUpload

# Add scripts/kaggle to path for drive_auth import
scripts_kaggle = Path(__file__).parent.parent / "kaggle"
sys.path.insert(0, str(scripts_kaggle))
import drive_auth

# Drive folder IDs
DRIVE_PUSH_FOLDER = drive_auth.DEFAULT_FOLDERS["push"]


def setup_drive(creds_path=None):
    """Authenticate to Drive using OAuth."""
    service, _ = drive_auth.setup_drive(creds_path=creds_path, verify=False)
    return service


def clear_push_folder(service):
    """Delete all files in push folder."""
    print("Clearing push folder...")
    results = service.files().list(
        q=f"'{DRIVE_PUSH_FOLDER}' in parents and trashed=false",
        fields="files(id, name)"
    ).execute()

    files = results.get("files", [])
    for f in files:
        service.files().delete(fileId=f["id"]).execute()
        print(f"  ✗ Deleted: {f['name']}")


def create_folder(service, name, parent_id):
    """Create folder in Drive."""
    meta = {
        "name": name,
        "mimeType": "application/vnd.google-apps.folder",
        "parents": [parent_id]
    }
    folder = service.files().create(body=meta, fields="id").execute()
    return folder["id"]


def upload_file(service, local_path, parent_id):
    """Upload file to Drive."""
    name = local_path.name
    size = local_path.stat().st_size
    size_str = f"{size:,} bytes" if size < 1024 else f"{size/1024:.1f} KB" if size < 1024*1024 else f"{size/(1024*1024):.1f} MB"

    meta = {
        "name": name,
        "parents": [parent_id]
    }
    media = MediaFileUpload(str(local_path))
    service.files().create(body=meta, media_body=media, fields="id").execute()
    print(f"    ↑ {name:40s} ({size_str})")


def build_push_directory(receptors, ligands, boxes, output_dir="kaggle_push", creds_path=None):
    """
    Upload all files needed for docking to Drive push folder.

    Args:
        receptors: List of receptor filenames (e.g., ["6X3U_apo.pdbqt"])
        ligands: List of ligand filenames (e.g., ["flumazenil.sdf"])
        boxes: List of box dicts with center/size (one per receptor)
        output_dir: Unused (kept for compatibility)
        creds_path: Unused (kept for compatibility, uses config.json OAuth)
    """
    # root is now the project root (up two levels)
    root = Path(__file__).parent.parent.parent
    service = setup_drive(creds_path)

    # Clear existing files
    clear_push_folder(service)

    print(f"\n=== Uploading to Drive push folder ===")

    # Create subdirectories
    receptors_folder_id = create_folder(service, "receptors", DRIVE_PUSH_FOLDER)
    ligands_folder_id = create_folder(service, "ligands", DRIVE_PUSH_FOLDER)
    print("  ✓ Created receptors/ and ligands/ folders")
    print()

    # Upload receptor files
    print("  Receptors:")
    receptor_dir = root / "data" / "receptors" / "prepped"
    for rec in receptors:
        src = receptor_dir / rec
        if not src.exists():
            raise FileNotFoundError(f"Receptor not found: {src}")
        upload_file(service, src, receptors_folder_id)

    print()
    # Upload ligand files
    print("  Ligands:")
    ligand_dir = root / "data" / "ligands"
    for lig in ligands:
        src = ligand_dir / lig
        if not src.exists():
            raise FileNotFoundError(f"Ligand not found: {src}")
        upload_file(service, src, ligands_folder_id)

    print()

    # Build and upload work queue with metadata
    jobs = []
    for rec, box in zip(receptors, boxes):
        for lig in ligands:
            jobs.append({
                "receptor": rec,
                "ligand": lig,
                "cx": box["center"][0],
                "cy": box["center"][1],
                "cz": box["center"][2],
                "sx": box["size"][0],
                "sy": box["size"][1],
                "sz": box["size"][2]
            })

    work_queue_data = {
        "task_type": "gnina_docking",
        "gnina_profile": "validation",
        "max_workers": 2,
        "jobs": jobs
    }

    print("  Work Queue:")
    print(f"    Task Type: {work_queue_data['task_type']}")
    print(f"    GNINA Profile: {work_queue_data['gnina_profile']}")
    print(f"    Max Workers: {work_queue_data['max_workers']}")
    print(f"    Jobs: {len(jobs)}")
    if jobs:
        print(f"    First Job: {jobs[0]['receptor']} + {jobs[0]['ligand']}")
    print()

    work_queue_file = root / "work_queue.json"
    with open(work_queue_file, 'w') as f:
        json.dump(work_queue_data, f, indent=2)

    upload_file(service, work_queue_file, DRIVE_PUSH_FOLDER)
    work_queue_file.unlink()  # Delete temp file

    # Upload scripts to scripts/ folder
    print()
    print("  Scripts:")
    scripts_folder_id = create_folder(service, "scripts", DRIVE_PUSH_FOLDER)

    kaggle_scripts = root / "scripts" / "kaggle"
    kaggle_script_files = [
        "setup_env.py",
        "drive_auth.py",
        "drive_io.py",
        "worker.py"
    ]

    for script in kaggle_script_files:
        src = kaggle_scripts / script
        if src.exists():
            upload_file(service, src, scripts_folder_id)

    # Upload gnina_runner.py from pipeline to scripts/ (flat structure)
    pipeline_scripts = root / "scripts" / "pipeline"
    gnina_runner = pipeline_scripts / "gnina_runner.py"
    if gnina_runner.exists():
        upload_file(service, gnina_runner, scripts_folder_id)

    print(f"\n{'='*60}")
    print(f"✓ Push folder ready on Drive")
    print(f"{'='*60}")
    print(f"  Jobs: {len(jobs)}")
    print(f"  Receptors: {len(receptors)}")
    print(f"  Ligands: {len(ligands)}")
    print(f"  Scripts: {len(kaggle_script_files) + 1}")
    print(f"\nNext step:")
    print(f"  cd scripts/kaggle && kaggle kernels push")


def main():
    parser = argparse.ArgumentParser(description="Upload docking job to Drive push folder")
    parser.add_argument("--receptor", nargs="+", help="Receptor filename(s) from data/receptors/prepped/")
    parser.add_argument("--ligand", nargs="+", help="Ligand filename(s) from data/ligands/")
    parser.add_argument("--box", help="Box params: cx,cy,cz,sx,sy,sz (e.g., 0,0,0,20,20,20)")
    parser.add_argument("--config", type=Path, help="JSON config file with receptors/ligands/boxes")

    args = parser.parse_args()

    if args.config:
        # Load from config file
        with open(args.config) as f:
            config = json.load(f)
        receptors = config["receptors"]
        ligands = config["ligands"]
        boxes = config["boxes"]
    elif args.receptor and args.ligand and args.box:
        # Load from CLI args
        receptors = args.receptor
        ligands = args.ligand
        box_vals = [float(x) for x in args.box.split(",")]
        boxes = [{
            "center": [box_vals[0], box_vals[1], box_vals[2]],
            "size": [box_vals[3], box_vals[4], box_vals[5]]
        }] * len(receptors)
    else:
        parser.error("Provide either --config or (--receptor + --ligand + --box)")

    build_push_directory(receptors, ligands, boxes)


if __name__ == "__main__":
    main()
