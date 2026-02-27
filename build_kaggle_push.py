#!/usr/bin/env python3
"""Upload files to Drive push folder for Kaggle docking jobs.

Usage:
    python build_kaggle_push.py --receptor 6X3U_apo.pdbqt --ligand flumazenil.sdf --box 0,0,0,20,20,20
    python build_kaggle_push.py --config job_config.json
"""

import json
import shutil
import argparse
from pathlib import Path
from google.oauth2 import service_account
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload

# Drive folder IDs
DRIVE_PUSH_FOLDER = "1_rJoAgziylRUQhWqXB2deTU1YW5_hBXM"


def setup_drive(creds_path):
    """Authenticate to Drive."""
    creds = service_account.Credentials.from_service_account_file(
        creds_path, scopes=["https://www.googleapis.com/auth/drive"]
    )
    return build("drive", "v3", credentials=creds)


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
    meta = {
        "name": name,
        "parents": [parent_id]
    }
    media = MediaFileUpload(str(local_path))
    service.files().create(body=meta, media_body=media, fields="id").execute()
    print(f"  ✓ Uploaded: {name}")


def build_push_directory(receptors, ligands, boxes, output_dir="kaggle_push", creds_path="service-account-key.json"):
    """
    Upload all files needed for docking to Drive push folder.

    Args:
        receptors: List of receptor filenames (e.g., ["6X3U_apo.pdbqt"])
        ligands: List of ligand filenames (e.g., ["flumazenil.sdf"])
        boxes: List of box dicts with center/size (one per receptor)
        output_dir: Unused (kept for compatibility)
        creds_path: Path to Google Drive service account key
    """
    root = Path(__file__).parent
    service = setup_drive(root / creds_path)

    # Clear existing files
    clear_push_folder(service)

    print(f"\nUploading to Drive push folder")

    # Create subdirectories
    receptors_folder_id = create_folder(service, "receptors", DRIVE_PUSH_FOLDER)
    ligands_folder_id = create_folder(service, "ligands", DRIVE_PUSH_FOLDER)
    print("  ✓ Created receptors/ and ligands/ folders")

    # Upload receptor files
    receptor_dir = root / "data" / "receptors" / "prepped"
    for rec in receptors:
        src = receptor_dir / rec
        if not src.exists():
            raise FileNotFoundError(f"Receptor not found: {src}")
        upload_file(service, src, receptors_folder_id)

    # Upload ligand files
    ligand_dir = root / "data" / "ligands"
    for lig in ligands:
        src = ligand_dir / lig
        if not src.exists():
            raise FileNotFoundError(f"Ligand not found: {src}")
        upload_file(service, src, ligands_folder_id)

    # Build and upload work queue
    work_queue = []
    for rec, box in zip(receptors, boxes):
        for lig in ligands:
            work_queue.append({
                "receptor": rec,
                "ligand": lig,
                "cx": box["center"][0],
                "cy": box["center"][1],
                "cz": box["center"][2],
                "sx": box["size"][0],
                "sy": box["size"][1],
                "sz": box["size"][2]
            })

    work_queue_file = root / "work_queue.json"
    with open(work_queue_file, 'w') as f:
        json.dump(work_queue, f, indent=2)

    upload_file(service, work_queue_file, DRIVE_PUSH_FOLDER)
    work_queue_file.unlink()  # Delete temp file
    print(f"  ✓ Uploaded work_queue.json ({len(work_queue)} jobs)")

    # Upload Kaggle scripts
    kaggle_scripts = root / "scripts" / "kaggle"
    needed_scripts = [
        "setup_env.py",
        "drive_auth.py",
        "drive_io.py",
        "gnina_worker.py"
    ]

    for script in needed_scripts:
        src = kaggle_scripts / script
        if src.exists():
            upload_file(service, src, DRIVE_PUSH_FOLDER)

    print(f"\n✓ Push folder ready on Drive")
    print(f"  Jobs: {len(work_queue)}")
    print(f"  Receptors: {len(receptors)}")
    print(f"  Ligands: {len(ligands)}")
    print(f"\nNow push bootstrap to Kaggle:")
    print(f"  cd kaggle_push && kaggle kernels push -p . --accelerator NvidiaTeslaT4")


def main():
    parser = argparse.ArgumentParser(description="Upload docking job to Drive push folder")
    parser.add_argument("--receptor", nargs="+", help="Receptor filename(s) from data/receptors/prepped/")
    parser.add_argument("--ligand", nargs="+", help="Ligand filename(s) from data/ligands/")
    parser.add_argument("--box", help="Box params: cx,cy,cz,sx,sy,sz (e.g., 0,0,0,20,20,20)")
    parser.add_argument("--config", type=Path, help="JSON config file with receptors/ligands/boxes")
    parser.add_argument("--creds", default="service-account-key.json", help="Path to Drive credentials")

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

    build_push_directory(receptors, ligands, boxes, "", args.creds)


if __name__ == "__main__":
    main()
