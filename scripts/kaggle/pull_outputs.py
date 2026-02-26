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
        "--creds",
        type=str,
        default="service-account-key.json",
        help="Path to service account key (default: service-account-key.json)"
    )
    parser.add_argument(
        "--folder",
        type=str,
        choices=["outputs", "validation", "receptors", "ligands"],
        default="outputs",
        help="Which Drive folder to download (default: outputs)"
    )
    args = parser.parse_args()

    # Resolve project root
    project_root = Path(__file__).parent.parent.parent
    creds_path = project_root / args.creds
    output_dir = project_root / args.output_dir

    if not creds_path.exists():
        print(f"✗ ERROR: Credentials not found at {creds_path}")
        print()
        print("Make sure service-account-key.json is in the project root.")
        sys.exit(1)

    print("=" * 60)
    print(f"Pulling {args.folder} from Google Drive")
    print("=" * 60)
    print()

    # Authenticate
    print("Authenticating...")
    drive_service, _ = drive_auth.setup_drive(str(creds_path), verify=False)
    print("✓ Authenticated")
    print()

    # Get folder ID
    folder_id = drive_auth.DEFAULT_FOLDERS[args.folder]

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Download
    print(f"Downloading to: {output_dir}")
    print()
    drive_io.download_folder(drive_service, folder_id, str(output_dir))

    print()
    print("=" * 60)
    print("✓ Download complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
