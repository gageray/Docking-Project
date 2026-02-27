#!/usr/bin/env python3
"""Bootstrap script for Kaggle notebooks - handles setup, execution, and cleanup."""
import subprocess
import sys
import os
import json
import shutil
from datetime import datetime
from pathlib import Path

# Drive folder IDs
DRIVE_PUSH_FOLDER = "1_rJoAgziylRUQhWqXB2deTU1YW5_hBXM"
DRIVE_OUTPUTS_FOLDER = "1u3U0R069wtYTLSBn2dRRaSXF74JPsoHk"

def get_credentials_path():
    """Find credentials in Kaggle dataset or working directory."""
    if os.path.exists("/kaggle/input/gdrive-oauth/key.json"):
        return "/kaggle/input/gdrive-oauth/key.json"
    elif os.path.exists("/kaggle/working/key.json"):
        return "/kaggle/working/key.json"
    return None

def download_workspace():
    """Download push folder from Google Drive."""
    print("\n=== Downloading workspace from Drive ===")

    creds_path = get_credentials_path()
    if not creds_path:
        raise Exception("No Drive credentials found!")

    from google.oauth2 import service_account
    from googleapiclient.discovery import build
    from googleapiclient.http import MediaIoBaseDownload

    creds = service_account.Credentials.from_service_account_file(
        creds_path, scopes=["https://www.googleapis.com/auth/drive"]
    )
    service = build("drive", "v3", credentials=creds)

    # Download all files from push folder
    results = service.files().list(
        q=f"'{DRIVE_PUSH_FOLDER}' in parents and trashed=false",
        fields="files(id, name, mimeType)"
    ).execute()

    files = results.get("files", [])
    for file in files:
        if file["mimeType"] == "application/vnd.google-apps.folder":
            download_folder(service, file["id"], file["name"], "/kaggle/working")
        else:
            download_file(service, file["id"], file["name"], "/kaggle/working")

    print(f"  ✓ Downloaded {len(files)} items from push folder")

def download_file(service, file_id, file_name, dest_dir):
    """Download a single file from Drive."""
    from googleapiclient.http import MediaIoBaseDownload
    request = service.files().get_media(fileId=file_id)
    dest_path = os.path.join(dest_dir, file_name)
    with open(dest_path, "wb") as f:
        downloader = MediaIoBaseDownload(f, request)
        done = False
        while not done:
            _, done = downloader.next_chunk()
    print(f"    ↓ {file_name}")

def download_folder(service, folder_id, folder_name, parent_dir):
    """Download a folder and its contents recursively."""
    dest_dir = os.path.join(parent_dir, folder_name)
    os.makedirs(dest_dir, exist_ok=True)
    results = service.files().list(
        q=f"'{folder_id}' in parents and trashed=false",
        fields="files(id, name, mimeType)"
    ).execute()
    files = results.get("files", [])
    for file in files:
        if file["mimeType"] == "application/vnd.google-apps.folder":
            download_folder(service, file["id"], file["name"], dest_dir)
        else:
            download_file(service, file["id"], file["name"], dest_dir)

def setup_environment():
    """Install dependencies and verify environment."""
    print("\n=== Setting up environment ===")
    sys.path.insert(0, "/kaggle/working")

    import setup_env
    setup_env.setup_kaggle_environment(verify=True)

def run_worker():
    """Execute the docking worker."""
    print("\n=== Running worker ===")
    os.chdir("/kaggle/working")

    import gnina_worker
    gnina_worker.main()

def generate_run_id():
    """Generate run ID in format YY-MM-DD-XX (hex increment)."""
    now = datetime.now()
    date_prefix = now.strftime("%y-%m-%d")

    creds_path = get_credentials_path()
    if not creds_path:
        return f"{date_prefix}-00"

    # Check Drive for existing runs today
    try:
        from google.oauth2 import service_account
        from googleapiclient.discovery import build

        creds = service_account.Credentials.from_service_account_file(
            creds_path, scopes=["https://www.googleapis.com/auth/drive"]
        )
        service = build("drive", "v3", credentials=creds)

        results = service.files().list(
            q=f"'{DRIVE_OUTPUTS_FOLDER}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
            fields="files(id, name)"
        ).execute()

        files = results.get("files", [])

        # Find max hex suffix for today
        max_suffix = -1
        for f in files:
            if f['name'].startswith(date_prefix):
                try:
                    suffix = f['name'].split('-')[-1]
                    max_suffix = max(max_suffix, int(suffix, 16))
                except:
                    pass

        next_suffix = max_suffix + 1
        return f"{date_prefix}-{next_suffix:02x}"

    except Exception as e:
        print(f"  ⚠ Could not check Drive for run ID, using 00: {e}")
        return f"{date_prefix}-00"

def create_run_metadata(run_id):
    """Create run metadata file."""
    metadata = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "worker": "gnina_worker.py"
    }

    # Copy work queue if exists
    if os.path.exists("/kaggle/working/work_queue.json"):
        with open("/kaggle/working/work_queue.json") as f:
            metadata["work_queue"] = json.load(f)

    # Write metadata
    with open("/kaggle/working/outputs/run_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"  ✓ Created run metadata: {run_id}")

def upload_outputs(run_id, zip_output=False):
    """Upload outputs folder to Drive."""
    print("\n=== Uploading outputs to Drive ===")

    creds_path = get_credentials_path()
    if not creds_path:
        print("  ⚠ No credentials found, skipping upload")
        return

    try:
        from google.oauth2 import service_account
        from googleapiclient.discovery import build
        from googleapiclient.http import MediaFileUpload

        creds = service_account.Credentials.from_service_account_file(
            creds_path, scopes=["https://www.googleapis.com/auth/drive"]
        )
        service = build("drive", "v3", credentials=creds)
        outputs_folder_id = DRIVE_OUTPUTS_FOLDER

        # Create run folder on Drive
        folder_meta = {
            "name": run_id,
            "mimeType": "application/vnd.google-apps.folder",
            "parents": [outputs_folder_id]
        }
        run_folder = service.files().create(body=folder_meta, fields="id").execute()
        run_folder_id = run_folder["id"]
        print(f"  ✓ Created Drive folder: {run_id}")

        # Upload all files from outputs/
        output_dir = Path("/kaggle/working/outputs")
        for file_path in output_dir.iterdir():
            if file_path.is_file():
                file_meta = {
                    "name": file_path.name,
                    "parents": [run_folder_id]
                }
                media = MediaFileUpload(str(file_path))
                service.files().create(
                    body=file_meta,
                    media_body=media,
                    fields="id"
                ).execute()
                print(f"  ✓ Uploaded: {file_path.name}")

        print(f"  ✓ All outputs uploaded to Drive: {run_id}")

    except Exception as e:
        print(f"  ✗ Upload failed: {e}")

def main():
    print("=" * 60)
    print("KAGGLE DOCKING BOOTSTRAP")
    print("=" * 60)

    try:
        # Startup
        download_workspace()
        setup_environment()

        # Generate run ID
        run_id = generate_run_id()
        print(f"\n=== Run ID: {run_id} ===")

        # Execute worker
        run_worker()

        # Shutdown
        create_run_metadata(run_id)
        upload_outputs(run_id)

        print("\n" + "=" * 60)
        print("BOOTSTRAP COMPLETE")
        print("=" * 60)

    except Exception as e:
        print(f"\n✗ Bootstrap failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
