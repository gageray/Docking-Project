# Kaggle Docking Utilities

Modular Python scripts for running molecular docking workflows on Kaggle notebooks with GPU acceleration and Google Drive persistence.

## Overview

This package provides reusable functions for:
- **Environment Setup** - Install dependencies (GNINA, RDKit, OpenMM, etc.)
- **Google Drive I/O** - Upload/download files with OAuth authentication
- **Repository Cloning** - Clone your docking pipeline into Kaggle
- **Directory Management** - Create standardized folder structures
- **Hardware Testing** - Verify GPU isolation for parallel docking

## Installation in Kaggle Notebook

```python
# Clone your docking repository
from scripts.kaggle import setup_repo

repo = setup_repo(
    repo_url="https://github.com/yourusername/docking-pipeline.git",
    dest_path="/kaggle/working/docking-pipeline",
    branch="main"
)
```

## Quick Start

### 1. Environment Setup

```python
from scripts.kaggle import setup_kaggle_environment

# Install all dependencies and verify
status = setup_kaggle_environment(verify=True)

# Check results
if status['gnina']['available']:
    print("✓ GNINA ready")
if all(pkg['available'] for pkg in status['packages'].values()):
    print("✓ All packages installed")
```

### 2. Google Drive Authentication

Upload your service account credentials as a Kaggle dataset first, then:

```python
from scripts.kaggle import setup_drive, DEFAULT_FOLDERS

# Authenticate and verify folders
drive_service, status = setup_drive(
    creds_path="/kaggle/input/your-dataset/service-account-key.json",
    verify=True
)

# Check access
if all(s['accessible'] for s in status.values()):
    print("✓ All Drive folders accessible")
```

### 3. Directory Setup

```python
from scripts.kaggle import setup_kaggle_dirs

# Create standard directory structure
dirs = setup_kaggle_dirs(
    base_path="/kaggle/working/docking",
    show_tree=True
)

print(dirs['receptors'])  # /kaggle/working/docking/receptors
print(dirs['ligands'])    # /kaggle/working/docking/ligands
```

### 4. Download Files from Drive

```python
from scripts.kaggle import download_folder, DEFAULT_FOLDERS

# Download all receptors
receptor_files = download_folder(
    drive_service,
    folder_id=DEFAULT_FOLDERS['receptors'],
    local_dir=dirs['receptors']
)

print(f"Downloaded {len(receptor_files)} receptor files")
```

### 5. Run Your Pipeline

```python
# Your docking code here
# e.g., from your cloned repo:
# import sys
# sys.path.insert(0, str(repo['path']))
# from scripts.receptor import receptor_prep
# ...
```

### 6. Upload Results to Drive

```python
from scripts.kaggle import upload_folder, DEFAULT_FOLDERS

# Upload all output files
uploaded_ids = upload_folder(
    drive_service,
    local_dir=dirs['outputs'],
    folder_id=DEFAULT_FOLDERS['outputs'],
    pattern="*.sdf"  # Only SDF files
)

print(f"Uploaded {len(uploaded_ids)} result files")
```

### 7. Hardware Verification (Optional)

```python
from scripts.kaggle import verify_hardware_setup

# Check GPU availability and test parallel isolation
hw_status = verify_hardware_setup(test_parallel=True)

print(f"GPUs available: {hw_status['gpus']['count']}")
print(f"Optimal workers: {hw_status['optimal_workers']}")
```

## Complete Notebook Example

```python
# ===== CELL 1: Setup Environment =====
from scripts.kaggle import setup_kaggle_environment
status = setup_kaggle_environment(verify=True)

# ===== CELL 2: Authenticate Drive =====
from scripts.kaggle import setup_drive
drive_service, _ = setup_drive(
    creds_path="/kaggle/input/creds/service-account-key.json"
)

# ===== CELL 3: Clone Repository =====
from scripts.kaggle import setup_repo
repo = setup_repo(
    repo_url="https://github.com/user/docking-pipeline.git",
    dest_path="/kaggle/working/docking-pipeline"
)

# ===== CELL 4: Create Directories =====
from scripts.kaggle import setup_kaggle_dirs
dirs = setup_kaggle_dirs()

# ===== CELL 5: Download Input Files =====
from scripts.kaggle import download_folder, DEFAULT_FOLDERS
download_folder(drive_service, DEFAULT_FOLDERS['receptors'], dirs['receptors'])
download_folder(drive_service, DEFAULT_FOLDERS['ligands'], dirs['ligands'])

# ===== CELL 6: Run Docking =====
# Your pipeline code here

# ===== CELL 7: Upload Results =====
from scripts.kaggle import upload_folder
upload_folder(drive_service, dirs['outputs'], DEFAULT_FOLDERS['outputs'])
```

## Module Reference

### `drive_auth.py`
- `authenticate_drive(creds_path)` - Create authenticated Drive service
- `verify_drive_folders(drive_service, folder_ids)` - Check folder access
- `setup_drive(creds_path, verify=True)` - Complete auth + verification

### `drive_io.py`
- `upload_to_drive(service, local_path, folder_id)` - Upload single file
- `download_from_drive(service, file_id, local_path)` - Download single file
- `list_drive_folder(service, folder_id)` - List folder contents
- `upload_folder(service, local_dir, folder_id, pattern)` - Upload multiple files
- `download_folder(service, folder_id, local_dir)` - Download multiple files

### `setup_env.py`
- `install_apt_packages(packages)` - Install system packages
- `install_pip_packages(packages)` - Install Python packages
- `install_gnina(install_path)` - Download GNINA binary
- `verify_environment()` - Check GPU, GNINA, packages
- `setup_kaggle_environment(verify=True)` - Complete setup

### `setup_dirs.py`
- `create_docking_dirs(base_path, subdirs)` - Create directory tree
- `show_directory_tree(base_path)` - Display structure
- `setup_kaggle_dirs(base_path, show_tree=True)` - Complete setup

### `clone_repo.py`
- `clone_repo(repo_url, dest_path, branch, token)` - Clone repository
- `pull_latest(repo_path, branch)` - Update existing repo
- `get_repo_info(repo_path)` - Get branch/commit info
- `setup_repo(repo_url, ...)` - Clone or update

### `hardware_check.py`
- `check_gpu_info()` - Query GPU count and specs
- `test_gpu_isolation(gpu_id)` - Test CUDA_VISIBLE_DEVICES
- `test_parallel_gpu_isolation(num_workers)` - Test parallel isolation
- `create_gpu_queue(gpu_ids)` - Create thread-safe GPU queue
- `get_optimal_workers()` - Recommend parallel worker count
- `verify_hardware_setup(test_parallel=True)` - Complete verification

## Drive Folder IDs

Default folder IDs are defined in `drive_auth.DEFAULT_FOLDERS`:

```python
{
    "root": "1NQ_5DucfsWMBRnZM6fSN82isem9Ha72R",
    "receptors": "1H3zFxbwk_GsTu65vmw8ZIsC6ykQRsNFI",
    "ligands": "1Y-9eDS9PgArkAllh3VrApTebqdWYLATH",
    "outputs": "1_NG8tk4Cz5bhQXvMctHhIDUIWUzChsUn",
    "validation": "1m9YVLgDzjX8DCCB-Zojej_inIxHd9vEJ"
}
```

You can override these when calling functions.

## Notes

- All functions include detailed docstrings with examples
- Functions handle errors gracefully and print status messages
- Service account credentials must be uploaded as a Kaggle dataset
- For private repos, pass `token` and `username` to `clone_repo()`
- GPU isolation tests require PyTorch (included in Kaggle by default)

## Future Kaggle Projects

These utilities are generic and can be reused for any Kaggle project requiring:
- Repository cloning
- Google Drive persistence
- Environment setup
- Parallel GPU execution
