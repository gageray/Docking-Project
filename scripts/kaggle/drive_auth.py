"""Google Drive authentication for Kaggle notebooks.

Handles OAuth2 authentication using refresh tokens.
"""

import os
import json
import requests
from google.oauth2.credentials import Credentials
from googleapiclient.discovery import build


# Default folder IDs for the docking project structure
DEFAULT_FOLDERS = {
    "root": "1NQ_5DucfsWMBRnZM6fSN82isem9Ha72R",       # docking/
    "push": "1_rJoAgziylRUQhWqXB2deTU1YW5_hBXM",       # docking/push/ (staged data for Kaggle)
    "receptors": "1H3zFxbwk_GsTu65vmw8ZIsC6ykQRsNFI",  # docking/receptors/
    "ligands": "1Y-9eDS9PgArkAllh3VrApTebqdWYLATH",    # docking/ligands/
    "outputs": "1_NG8tk4Cz5bhQXvMctHhIDUIWUzChsUn",    # docking/outputs/
    "validation": "1m9YVLgDzjX8DCCB-Zojej_inIxHd9vEJ"  # docking/validation/
}


def get_oauth_credentials(creds_path=None):
    """Load OAuth credentials from drive_oauth.json or config.json.

    Returns:
        dict: {client_id, client_secret, refresh_token}
    """
    # 1. Try explicit creds_path or default drive_oauth.json in root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    if creds_path is None:
        creds_path = os.path.join(project_root, 'drive_oauth.json')
    
    if os.path.exists(creds_path):
        with open(creds_path, 'r') as f:
            content = f.read().strip()
        
        # Handle the "drive_oauth": { ... } format without outer braces
        if content.startswith('"drive_oauth"'):
            try:
                data = json.loads(f"{{{content}}}")
                return data['drive_oauth']
            except json.JSONDecodeError:
                pass
        
        try:
            data = json.loads(content)
            return data.get('drive_oauth', data)
        except json.JSONDecodeError:
            pass

    # 3. Fall back to config.json
    config_path = os.path.join(project_root, 'config.json')
    if os.path.exists(config_path):
        with open(config_path) as f:
            config = json.load(f)
        if 'drive_oauth' in config:
            return config['drive_oauth']
            
    raise RuntimeError(f"Could not find valid credentials at {creds_path} or in config.json")


def refresh_access_token(client_id, client_secret, refresh_token):
    """Exchange refresh token for a new access token.

    Args:
        client_id (str): OAuth client ID
        client_secret (str): OAuth client secret
        refresh_token (str): Refresh token

    Returns:
        str: New access token
    """
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
    return response.json()['access_token']


def authenticate_drive(creds_path=None, scopes=None):
    """Authenticate to Google Drive using OAuth2.

    Args:
        creds_path (str, optional): Path to OAuth credentials.
        scopes (list, optional): OAuth scopes. Defaults to Drive access.

    Returns:
        Resource: Authenticated Drive service object

    Example:
        >>> drive_service = authenticate_drive()
    """
    if scopes is None:
        scopes = ["https://www.googleapis.com/auth/drive"]

    oauth_config = get_oauth_credentials(creds_path)
    access_token = refresh_access_token(
        oauth_config['client_id'],
        oauth_config['client_secret'],
        oauth_config['refresh_token']
    )

    creds = Credentials(token=access_token)
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


def setup_drive(creds_path=None, verify=True, folder_ids=None):
    """Complete Drive setup: authenticate and optionally verify folders.

    Args:
        creds_path (str, optional): Path to OAuth credentials.
        verify (bool): Whether to verify folder access. Default True.
        folder_ids (dict, optional): Custom folder IDs to verify

    Returns:
        tuple: (drive_service, verification_results or None)

    Example:
        >>> service, status = setup_drive()
        >>> if all(s['accessible'] for s in status.values()):
        ...     print("All folders ready!")
    """
    drive_service = authenticate_drive(creds_path)

    if verify:
        verification = verify_drive_folders(drive_service, folder_ids)
        return drive_service, verification

    return drive_service, None
