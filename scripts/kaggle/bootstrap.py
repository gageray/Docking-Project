#!/usr/bin/env python3
"""Bootstrap script for Kaggle notebooks - handles setup, execution, and cleanup.

IMPORTANT: To push this kernel with T4x2 GPUs instead of P100, use:
    kaggle kernels push --accelerator NvidiaTeslaT4

Available accelerators: NvidiaTeslaP100 (default), NvidiaTeslaT4, NvidiaTeslaT4Highmem, TpuV6E8

OAuth credentials are loaded from private dataset: ineptrobot/drive-oauth
Mounted at: /kaggle/input/drive-oauth/oauth_credentials.json
"""
import subprocess
import sys
import os
import json
import shutil
from datetime import datetime
from pathlib import Path
import requests
from google.oauth2.credentials import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaIoBaseDownload, MediaFileUpload

# Drive folder IDs
DEFAULT_FOLDERS = {
    "push": "1_rJoAgziylRUQhWqXB2deTU1YW5_hBXM",
    "outputs": "1_NG8tk4Cz5bhQXvMctHhIDUIWUzChsUn",
}

def authenticate_drive():
    """Authenticate to Google Drive using private dataset credentials.
    
    Self-contained in bootstrap.py to allow initial workspace download before 
    the rest of the scripts are available.
    """
    # Path is /kaggle/input/[dataset-slug]/[filename]
    # The dataset slug is drive-oauth, but the file is drive_oauth.json
    SECRET_PATH = "/kaggle/input/drive-oauth/drive_oauth.json"

    if not os.path.exists(SECRET_PATH):
        # Fallback to check if the directory name also changed
        ALT_PATH = "/kaggle/input/drive_oauth/drive_oauth.json"
        if os.path.exists(ALT_PATH):
            SECRET_PATH = ALT_PATH
        else:
            raise RuntimeError(f"OAuth credentials not found. Checked {SECRET_PATH} and {ALT_PATH}. Ensure 'ineptrobot/drive-oauth' is in dataset_sources.")

    with open(SECRET_PATH, 'r') as f:
        content = f.read().strip()
    
    try:
        # Handle the "drive_oauth": { ... } format without outer braces
        if content.startswith('"drive_oauth"'):
            oauth_data = json.loads(f"{{{content}}}")["drive_oauth"]
        else:
            oauth_data = json.loads(content)
            if "drive_oauth" in oauth_data:
                oauth_data = oauth_data["drive_oauth"]
    except Exception as e:
        raise RuntimeError(f"Could not parse JSON from {SECRET_PATH}: {e}")

    client_id = oauth_data['client_id']
    client_secret = oauth_data['client_secret']
    refresh_token = oauth_data['refresh_token']

    # Get access token
    response = requests.post(
        'https://oauth2.googleapis.com/token',
        data={
            'client_id': client_id,
            'client_secret': client_secret,
            'refresh_token': refresh_token,
            'grant_type': 'refresh_token'
        }
    )
    response.raise_for_status()
    access_token = response.json()['access_token']

    creds = Credentials(token=access_token)
    return build("drive", "v3", credentials=creds)

def download_workspace():
    """Download push folder from Google Drive."""
    print("\n=== Downloading workspace from Drive ===")

    service = authenticate_drive()

    # Download all files from push folder
    push_folder_id = DEFAULT_FOLDERS["push"]
    print(f"  Source: push folder (ID: {push_folder_id[:8]}...)")

    results = service.files().list(
        q=f"'{push_folder_id}' in parents and trashed=false",
        fields="files(id, name, mimeType)"
    ).execute()

    files = results.get("files", [])
    print(f"  Found {len(files)} items in push folder")
    print()

    for file in files:
        if file["mimeType"] == "application/vnd.google-apps.folder":
            print(f"  [DIR] {file['name']}/")
            download_folder(service, file["id"], file["name"], "/kaggle/working")
        else:
            download_file(service, file["id"], file["name"], "/kaggle/working")

    print()
    print(f"  SUCCESS: Downloaded {len(files)} items from push folder")

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

    # Get file size for logging
    size = os.path.getsize(dest_path)
    size_str = f"{size:,} bytes" if size < 1024 else f"{size/1024:.1f} KB" if size < 1024*1024 else f"{size/(1024*1024):.1f} MB"
    print(f"    DOWNLOAD: {file_name:40s} ({size_str})")

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
    
    # Scripts are inside the 'scripts' folder
    SCRIPTS_DIR = "/kaggle/working/scripts"
    if os.path.exists(SCRIPTS_DIR):
        sys.path.insert(0, SCRIPTS_DIR)
    sys.path.insert(0, "/kaggle/working")

    import setup_env
    setup_env.setup_kaggle_environment(verify=True)

def print_work_queue():
    """Display work queue contents before running worker."""
    work_queue_path = "/kaggle/working/work_queue.json"
    if not os.path.exists(work_queue_path):
        print("  WARNING: work_queue.json not found!")
        return

    with open(work_queue_path) as f:
        data = json.load(f)

    print("\n=== Work Queue Loaded ===")
    print(f"  Task Type: {data.get('task_type', 'UNKNOWN')}")
    print(f"  GNINA Profile: {data.get('gnina_profile', 'UNKNOWN')}")
    print(f"  Max Workers: {data.get('max_workers', 'UNKNOWN')}")
    print(f"  Jobs: {len(data.get('jobs', []))}")

    # Show first job details
    jobs = data.get('jobs', [])
    if jobs:
        job = jobs[0]
        print(f"\n  First Job:")
        print(f"    Receptor: {job.get('receptor', 'N/A')}")
        print(f"    Ligand: {job.get('ligand', 'N/A')}")
        print(f"    Box Center: ({job.get('cx', 0)}, {job.get('cy', 0)}, {job.get('cz', 0)})")
        print(f"    Box Size: ({job.get('sx', 0)}, {job.get('sy', 0)}, {job.get('sz', 0)})")

def run_worker():
    """Execute the docking worker."""
    print_work_queue()

    print("\n=== Running worker ===")
    os.chdir("/kaggle/working")

    import worker
    worker.main()

def generate_run_id():
    """Generate run ID in format YY-MM-DD-XX (hex increment)."""
    now = datetime.now()
    date_prefix = now.strftime("%y-%m-%d")

    print("\n=== Generating Run ID ===")
    print(f"  Date prefix: {date_prefix}")

    # Check Drive for existing runs today
    try:
        service = authenticate_drive()

        outputs_folder_id = DEFAULT_FOLDERS["outputs"]
        results = service.files().list(
            q=f"'{outputs_folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
            fields="files(id, name)"
        ).execute()

        files = results.get("files", [])
        print(f"  Found {len(files)} existing run folders in Drive outputs/")

        # Find max hex suffix for today
        max_suffix = -1
        today_runs = []
        for f in files:
            if f['name'].startswith(date_prefix):
                today_runs.append(f['name'])
                try:
                    suffix = f['name'].split('-')[-1]
                    max_suffix = max(max_suffix, int(suffix, 16))
                except:
                    pass

        if today_runs:
            print(f"  Existing runs today: {', '.join(today_runs)}")
            print(f"  Max suffix: {max_suffix:02x}")

        next_suffix = max_suffix + 1
        run_id = f"{date_prefix}-{next_suffix:02x}"
        print(f"  Generated Run ID: {run_id}")
        return run_id

    except Exception as e:
        print(f"  WARNING: Could not check Drive for run ID, using 00: {e}")
        return f"{date_prefix}-00"

def create_run_metadata(run_id):
    """Create run metadata file."""
    metadata = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "worker": "worker.py"
    }

    # Copy work queue if exists
    if os.path.exists("/kaggle/working/work_queue.json"):
        with open("/kaggle/working/work_queue.json") as f:
            metadata["work_queue"] = json.load(f)

    # Write metadata
    with open("/kaggle/working/outputs/run_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"  SUCCESS: Created run metadata: {run_id}")

def upload_outputs(run_id, zip_output=False):
    """Upload outputs folder to Drive."""
    print("\n=== Uploading outputs to Drive ===")
    print(f"  Run ID: {run_id}")

    try:
        service = authenticate_drive()
        outputs_folder_id = DEFAULT_FOLDERS["outputs"]

        # Create run folder on Drive
        folder_meta = {
            "name": run_id,
            "mimeType": "application/vnd.google-apps.folder",
            "parents": [outputs_folder_id]
        }
        run_folder = service.files().create(body=folder_meta, fields="id").execute()
        run_folder_id = run_folder["id"]
        print(f"  SUCCESS: Created Drive folder: {run_id}")
        print()

        # Upload all files from outputs/
        output_dir = Path("/kaggle/working/outputs")
        files_uploaded = 0
        total_size = 0

        for file_path in output_dir.iterdir():
            if file_path.is_file():
                size = file_path.stat().st_size
                size_str = f"{size:,} bytes" if size < 1024 else f"{size/1024:.1f} KB" if size < 1024*1024 else f"{size/(1024*1024):.1f} MB"

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
                print(f"    UPLOAD: {file_path.name:40s} ({size_str})")

                files_uploaded += 1
                total_size += size

        print()
        total_size_str = f"{total_size:,} bytes" if total_size < 1024 else f"{total_size/1024:.1f} KB" if total_size < 1024*1024 else f"{total_size/(1024*1024):.1f} MB"
        print(f"  SUCCESS: Uploaded {files_uploaded} files ({total_size_str}) to Drive: {run_id}")

    except Exception as e:
        print(f"  FAILED: Upload failed: {e}")

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
        print(f"\nERROR: Bootstrap failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
