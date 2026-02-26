"""Google Drive I/O operations for Kaggle notebooks.

Upload/download files and list Drive folder contents.
"""

import os
from googleapiclient.http import MediaFileUpload, MediaIoBaseDownload


def upload_to_drive(drive_service, local_path, folder_id=None, overwrite=True):
    """Upload a local file to Google Drive.

    Args:
        drive_service: Authenticated Drive service object
        local_path (str): Path to local file to upload
        folder_id (str, optional): Parent folder ID. If None, uploads to root.
        overwrite (bool): If True and file exists, update it. Default True.

    Returns:
        str: File ID of uploaded/updated file

    Example:
        >>> file_id = upload_to_drive(service, "results.sdf", folder_id="abc123")
        >>> print(f"Uploaded: {file_id}")
    """
    name = os.path.basename(local_path)

    # Check for existing file
    q = f"name='{name}' and trashed=false"
    if folder_id:
        q += f" and '{folder_id}' in parents"

    existing = (
        drive_service.files()
        .list(q=q, fields="files(id)")
        .execute()
        .get("files", [])
    )

    media = MediaFileUpload(local_path)

    if existing and overwrite:
        # Update existing file
        updated = drive_service.files().update(
            fileId=existing[0]["id"], media_body=media
        ).execute()
        print(f"  ↑ updated {name} ({updated['id']})")
        return updated["id"]
    else:
        # Create new file
        meta = {"name": name}
        if folder_id:
            meta["parents"] = [folder_id]

        created = drive_service.files().create(
            body=meta, media_body=media, fields="id"
        ).execute()
        print(f"  ↑ created {name} ({created['id']})")
        return created["id"]


def download_from_drive(drive_service, file_id, local_path):
    """Download a file from Google Drive by ID.

    Args:
        drive_service: Authenticated Drive service object
        file_id (str): Drive file ID to download
        local_path (str): Local path to save file

    Example:
        >>> download_from_drive(service, "xyz789", "/kaggle/working/output.pdb")
    """
    request = drive_service.files().get_media(fileId=file_id)

    with open(local_path, "wb") as f:
        downloader = MediaIoBaseDownload(f, request)
        done = False
        while not done:
            _, done = downloader.next_chunk()

    print(f"  ↓ {local_path}")


def list_drive_folder(drive_service, folder_id, verbose=True):
    """List all files in a Drive folder.

    Args:
        drive_service: Authenticated Drive service object
        folder_id (str): Drive folder ID to list
        verbose (bool): Print file details. Default True.

    Returns:
        list: List of file dicts with 'id', 'name', 'mimeType', 'size'

    Example:
        >>> files = list_drive_folder(service, "abc123")
        >>> pdb_files = [f for f in files if f['name'].endswith('.pdb')]
    """
    results = drive_service.files().list(
        q=f"'{folder_id}' in parents and trashed=false",
        fields="files(id,name,mimeType,size)"
    ).execute()

    files = results.get("files", [])

    if verbose:
        for f in files:
            size = f.get("size", "folder")
            print(f"  {f['name']:40s}  {size:>10s}  {f['id']}")

    return files


def download_folder(drive_service, folder_id, local_dir):
    """Download all files from a Drive folder to a local directory.

    Args:
        drive_service: Authenticated Drive service object
        folder_id (str): Drive folder ID
        local_dir (str): Local directory path to save files

    Returns:
        list: Paths of downloaded files

    Example:
        >>> paths = download_folder(service, "xyz789", "/kaggle/working/receptors")
        >>> print(f"Downloaded {len(paths)} files")
    """
    os.makedirs(local_dir, exist_ok=True)

    files = list_drive_folder(drive_service, folder_id, verbose=False)
    downloaded = []

    for f in files:
        # Skip folders (they don't have a size attribute in the API)
        if "size" not in f:
            continue

        local_path = os.path.join(local_dir, f["name"])
        download_from_drive(drive_service, f["id"], local_path)
        downloaded.append(local_path)

    return downloaded


def upload_folder(drive_service, local_dir, folder_id, pattern="*"):
    """Upload all files from a local directory to a Drive folder.

    Args:
        drive_service: Authenticated Drive service object
        local_dir (str): Local directory containing files
        folder_id (str): Target Drive folder ID
        pattern (str): Glob pattern for files to upload. Default "*" (all).

    Returns:
        list: File IDs of uploaded files

    Example:
        >>> file_ids = upload_folder(service, "/kaggle/working/outputs",
        ...                           "abc123", pattern="*.sdf")
    """
    from pathlib import Path

    local_path = Path(local_dir)
    files = list(local_path.glob(pattern))

    uploaded = []
    for file_path in files:
        if file_path.is_file():
            file_id = upload_to_drive(drive_service, str(file_path), folder_id)
            uploaded.append(file_id)

    return uploaded
