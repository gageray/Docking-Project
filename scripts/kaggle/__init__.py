"""Kaggle utilities for molecular docking workflows.

Modular scripts for environment setup, Google Drive I/O, repository cloning,
and hardware management on Kaggle notebooks.

Modules:
    drive_auth: Google Drive authentication
    drive_io: Google Drive file upload/download
    setup_env: Install dependencies and verify environment
    setup_dirs: Create directory structure
    clone_repo: Git repository cloning
    hardware_check: GPU isolation and parallel execution testing
"""

from .drive_auth import authenticate_drive, verify_drive_folders, setup_drive, DEFAULT_FOLDERS
from .drive_io import (
    upload_to_drive,
    download_from_drive,
    list_drive_folder,
    download_folder,
    upload_folder
)
from .setup_env import (
    install_apt_packages,
    install_pip_packages,
    install_gnina,
    verify_environment,
    setup_kaggle_environment
)
from .setup_dirs import create_docking_dirs, show_directory_tree, setup_kaggle_dirs
from .clone_repo import clone_repo, pull_latest, get_repo_info, setup_repo
from .hardware_check import (
    check_gpu_info,
    test_gpu_isolation,
    test_parallel_gpu_isolation,
    create_gpu_queue,
    get_optimal_workers,
    verify_hardware_setup
)

__all__ = [
    # Drive auth
    "authenticate_drive",
    "verify_drive_folders",
    "setup_drive",
    "DEFAULT_FOLDERS",
    # Drive I/O
    "upload_to_drive",
    "download_from_drive",
    "list_drive_folder",
    "download_folder",
    "upload_folder",
    # Environment
    "install_apt_packages",
    "install_pip_packages",
    "install_gnina",
    "verify_environment",
    "setup_kaggle_environment",
    # Directories
    "create_docking_dirs",
    "show_directory_tree",
    "setup_kaggle_dirs",
    # Repository
    "clone_repo",
    "pull_latest",
    "get_repo_info",
    "setup_repo",
    # Hardware
    "check_gpu_info",
    "test_gpu_isolation",
    "test_parallel_gpu_isolation",
    "create_gpu_queue",
    "get_optimal_workers",
    "verify_hardware_setup",
]
