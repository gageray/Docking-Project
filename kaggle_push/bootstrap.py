#!/usr/bin/env python3
"""Bootstrap script for Kaggle notebooks - handles setup, execution, and cleanup."""
import subprocess
import sys
import os
import json
import shutil
from datetime import datetime
from pathlib import Path

REPO_URL = "https://github.com/gageray/Docking-Project.git"
SPARSE_PATH = "kaggle_push"

def run_cmd(cmd, check=True, capture=False):
    """Run shell command with error handling."""
    result = subprocess.run(
        cmd, shell=True, check=check,
        capture_output=capture, text=True
    )
    return result.stdout if capture else None

def clone_workspace():
    """Sparse clone just the kaggle_push folder."""
    print("\n=== Cloning workspace from GitHub ===")

    # Clone with sparse checkout
    run_cmd("git clone --depth 1 --filter=blob:none --sparse " + REPO_URL + " /tmp/repo")
    run_cmd("cd /tmp/repo && git sparse-checkout set " + SPARSE_PATH)

    # Copy kaggle_push contents to /kaggle/working
    src = Path(f"/tmp/repo/{SPARSE_PATH}")
    for item in src.iterdir():
        dest = Path("/kaggle/working") / item.name
        if item.is_dir():
            shutil.copytree(item, dest, dirs_exist_ok=True)
        else:
            shutil.copy2(item, dest)

    print(f"  ✓ Workspace cloned to /kaggle/working")

    # Cleanup
    shutil.rmtree("/tmp/repo")

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

    # Check if credentials exist for Drive check
    if not os.path.exists("/kaggle/working/key.json"):
        return f"{date_prefix}-00"

    # Check Drive for existing runs today
    try:
        sys.path.insert(0, "/kaggle/working")
        import drive_auth, drive_io

        service, _ = drive_auth.setup_drive("/kaggle/working/key.json", verify=False)
        files = drive_io.list_drive_folder(
            service,
            drive_auth.DEFAULT_FOLDERS["outputs"],
            verbose=False
        )

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

    if not os.path.exists("/kaggle/working/key.json"):
        print("  ⚠ No credentials found, skipping upload")
        return

    try:
        sys.path.insert(0, "/kaggle/working")
        import drive_auth
        from googleapiclient.http import MediaFileUpload

        service, _ = drive_auth.setup_drive("/kaggle/working/key.json", verify=False)
        outputs_folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]

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
        clone_workspace()
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
