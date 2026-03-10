#!/usr/bin/env python3
"""
Drive Utilities for Analysis Scripts
Helper functions for finding and downloading output folders from Google Drive.
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List

# Add kaggle scripts to path for Drive auth
sys.path.insert(0, str(Path(__file__).parent.parent / 'kaggle'))
import drive_auth
import drive_io


def find_most_recent_output(drive_service) -> Optional[Dict]:
    """
    Find the most recent output folder in Google Drive.

    Args:
        drive_service: Authenticated Drive service

    Returns:
        Dict with 'id', 'name', 'createdTime' or None if no outputs found
    """
    folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]

    runs = drive_service.files().list(
        q=f"'{folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
        fields="files(id, name, createdTime)",
        orderBy="createdTime desc"
    ).execute().get("files", [])

    return runs[0] if runs else None


def find_output_by_name(drive_service, folder_name: str) -> Optional[Dict]:
    """
    Find a specific output folder by name in Google Drive.

    Args:
        drive_service: Authenticated Drive service
        folder_name: Name of the output folder to find

    Returns:
        Dict with 'id', 'name', 'createdTime' or None if not found
    """
    folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]

    runs = drive_service.files().list(
        q=f"'{folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
        fields="files(id, name, createdTime)",
        orderBy="createdTime desc"
    ).execute().get("files", [])

    for run in runs:
        if run['name'] == folder_name:
            return run

    return None


def list_all_outputs(drive_service) -> List[Dict]:
    """
    List all output folders in Google Drive.

    Args:
        drive_service: Authenticated Drive service

    Returns:
        List of dicts with 'id', 'name', 'createdTime'
    """
    folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]

    runs = drive_service.files().list(
        q=f"'{folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and trashed=false",
        fields="files(id, name, createdTime)",
        orderBy="createdTime desc"
    ).execute().get("files", [])

    return runs


def search_outputs(drive_service, search_term: str) -> List[Dict]:
    """
    Search for output folders matching a term (e.g., date pattern).

    Args:
        drive_service: Authenticated Drive service
        search_term: String to search for in folder names

    Returns:
        List of matching folders
    """
    folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]

    runs = drive_service.files().list(
        q=f"'{folder_id}' in parents and mimeType='application/vnd.google-apps.folder' and name contains '{search_term}' and trashed=false",
        fields="files(id, name, createdTime)",
        orderBy="createdTime desc"
    ).execute().get("files", [])

    return runs


def download_output_folder(drive_service, folder_info: Dict, local_base_dir: str = "data/outputs") -> Path:
    """
    Download an output folder from Drive to local directory.

    Args:
        drive_service: Authenticated Drive service
        folder_info: Dict with 'id' and 'name' from find/search functions
        local_base_dir: Base directory for downloads

    Returns:
        Path to downloaded folder
    """
    project_root = Path(__file__).parent.parent.parent
    local_dir = project_root / local_base_dir / folder_info['name']

    if not local_dir.exists():
        local_dir.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {folder_info['name']} to {local_dir}...")
        drive_io.download_folder(drive_service, folder_info['id'], str(local_dir))
        print(f"Download complete!")
    else:
        print(f"Directory {local_dir} already exists. Skipping download.")

    return local_dir


def main():
    """CLI for browsing and downloading Drive outputs."""
    import argparse

    parser = argparse.ArgumentParser(description="Browse and download Google Drive outputs")
    parser.add_argument('--list', action='store_true', help='List all output folders')
    parser.add_argument('--search', type=str, help='Search for folders matching term')
    parser.add_argument('--download-recent', action='store_true', help='Download most recent output')
    parser.add_argument('--download-name', type=str, help='Download specific folder by name')
    parser.add_argument('--local-dir', type=str, default='data/outputs', help='Local directory for downloads')

    args = parser.parse_args()

    # Authenticate
    print("Authenticating to Google Drive...")
    drive_service, _ = drive_auth.setup_drive(verify=False)

    if args.list:
        print("\nAvailable output folders:")
        print("=" * 80)
        runs = list_all_outputs(drive_service)
        for i, run in enumerate(runs, 1):
            print(f"{i:3d}. {run['name']:40s} (created: {run['createdTime'][:10]})")
        print(f"\nTotal: {len(runs)} folders")

    elif args.search:
        print(f"\nSearching for: {args.search}")
        print("=" * 80)
        runs = search_outputs(drive_service, args.search)
        for i, run in enumerate(runs, 1):
            print(f"{i:3d}. {run['name']:40s} (created: {run['createdTime'][:10]})")
        print(f"\nFound: {len(runs)} folders")

    elif args.download_recent:
        run = find_most_recent_output(drive_service)
        if run:
            print(f"\nMost recent output: {run['name']} (created: {run['createdTime'][:10]})")
            local_dir = download_output_folder(drive_service, run, args.local_dir)
            print(f"\nDownloaded to: {local_dir}")
        else:
            print("No output folders found!")

    elif args.download_name:
        run = find_output_by_name(drive_service, args.download_name)
        if run:
            print(f"\nFound: {run['name']} (created: {run['createdTime'][:10]})")
            local_dir = download_output_folder(drive_service, run, args.local_dir)
            print(f"\nDownloaded to: {local_dir}")
        else:
            print(f"Folder '{args.download_name}' not found!")

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
