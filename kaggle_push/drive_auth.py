"""Google Drive authentication for Kaggle notebooks.

Handles OAuth setup using service account credentials from Kaggle datasets.
"""

from google.oauth2 import service_account
from googleapiclient.discovery import build


# Default folder IDs for the docking project structure
DEFAULT_FOLDERS = {
    "root": "1NQ_5DucfsWMBRnZM6fSN82isem9Ha72R",       # docking/
    "receptors": "1H3zFxbwk_GsTu65vmw8ZIsC6ykQRsNFI",  # docking/receptors/
    "ligands": "1Y-9eDS9PgArkAllh3VrApTebqdWYLATH",    # docking/ligands/
    "outputs": "1_NG8tk4Cz5bhQXvMctHhIDUIWUzChsUn",    # docking/outputs/
    "validation": "1m9YVLgDzjX8DCCB-Zojej_inIxHd9vEJ"  # docking/validation/
}


def authenticate_drive(creds_path, scopes=None):
    """Authenticate to Google Drive using service account credentials.

    Args:
        creds_path (str): Path to service account JSON key file
        scopes (list, optional): OAuth scopes. Defaults to Drive access.

    Returns:
        Resource: Authenticated Drive service object

    Example:
        >>> drive_service = authenticate_drive("/kaggle/input/my-creds/key.json")
    """
    if scopes is None:
        scopes = ["https://www.googleapis.com/auth/drive"]

    creds = service_account.Credentials.from_service_account_file(
        creds_path, scopes=scopes
    )
    drive_service = build("drive", "v3", credentials=creds)

    return drive_service


def verify_drive_folders(drive_service, folder_ids=None):
    """Verify that all Drive folders are accessible.

    Args:
        drive_service: Authenticated Drive service object
        folder_ids (dict, optional): Dict of {name: folder_id}. Defaults to project folders.

    Returns:
        dict: {folder_name: status_dict} with 'accessible' bool and 'name' or 'error'

    Example:
        >>> status = verify_drive_folders(drive_service)
        >>> if status['receptors']['accessible']:
        ...     print(f"Receptors folder: {status['receptors']['name']}")
    """
    if folder_ids is None:
        folder_ids = DEFAULT_FOLDERS

    results = {}

    for name, fid in folder_ids.items():
        try:
            f = drive_service.files().get(fileId=fid, fields="name").execute()
            results[name] = {"accessible": True, "name": f["name"]}
            print(f"  OK  {name:12s} -> {f['name']}")
        except Exception as e:
            results[name] = {"accessible": False, "error": str(e)}
            print(f"  FAIL {name:12s} -> {e}")

    return results


def setup_drive(creds_path, verify=True, folder_ids=None):
    """Complete Drive setup: authenticate and optionally verify folders.

    Args:
        creds_path (str): Path to service account JSON key
        verify (bool): Whether to verify folder access. Default True.
        folder_ids (dict, optional): Custom folder IDs to verify

    Returns:
        tuple: (drive_service, verification_results or None)

    Example:
        >>> service, status = setup_drive("/kaggle/input/creds/key.json")
        >>> if all(s['accessible'] for s in status.values()):
        ...     print("All folders ready!")
    """
    drive_service = authenticate_drive(creds_path)

    if verify:
        verification = verify_drive_folders(drive_service, folder_ids)
        return drive_service, verification

    return drive_service, None
