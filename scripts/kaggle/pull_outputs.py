#!/usr/bin/env python3
"""Pull docking outputs from Google Drive to local machine.

Useful for downloading results after Kaggle runs without using the browser.
"""

import sys
import os
from pathlib import Path

# Add kaggle scripts to path
sys.path.insert(0, str(Path(__file__).parent))

import drive_auth
import drive_io


def list_available_runs(service, folder_id):
    """List all run folders in outputs/ folder."""
    results = service.files().list(
        q=f"'{folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
        fields="files(id, name, createdTime)",
        orderBy="createdTime desc"
    ).execute()

    folders = results.get("files", [])
    return folders


def main():
    """Download outputs from Drive."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Download docking outputs from Google Drive"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data/outputs",
        help="Local directory to save outputs (default: data/outputs)"
    )
    parser.add_argument(
        "--folder",
        type=str,
        choices=["outputs", "validation", "receptors", "ligands"],
        default="outputs",
        help="Which Drive folder to download (default: outputs)"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available run folders without downloading"
    )
    args = parser.parse_args()

    # Resolve project root
    project_root = Path(__file__).parent.parent.parent
    output_dir = project_root / args.output_dir

    print("=" * 60)
    print(f"Pulling {args.folder} from Google Drive")
    print("=" * 60)
    print()

    # Authenticate using OAuth (same as rest of pipeline)
    print("Authenticating via OAuth...")
    drive_service, _ = drive_auth.setup_drive(verify=False)
    print("✓ Authenticated")
    print()

    # Get folder ID
    folder_id = drive_auth.DEFAULT_FOLDERS[args.folder]
    print(f"Target folder ID: {folder_id[:8]}...")
    print()

    # If listing mode, show available runs and exit
    if args.list:
        print("Available run folders:")
        print()
        runs = list_available_runs(drive_service, folder_id)
        if not runs:
            print("  No run folders found")
        else:
            for run in runs:
                print(f"  {run['name']:20s}  (created: {run['createdTime'][:10]})")
        print()
        print(f"Total: {len(runs)} runs")
        return

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download
    print(f"Downloading to: {output_dir}")
    print()

    # Show available runs
    runs = list_available_runs(drive_service, folder_id)
    if runs:
        print(f"Found {len(runs)} run folders:")
        for run in runs[:5]:  # Show first 5
            print(f"  - {run['name']}")
        if len(runs) > 5:
            print(f"  ... and {len(runs) - 5} more")
        print()

    drive_io.download_folder(drive_service, folder_id, str(output_dir))

    print()
    print("=" * 60)
    print("✓ Download complete!")
    print("=" * 60)
    print(f"Files saved to: {output_dir}")


if __name__ == "__main__":
    main()
