"""Directory structure setup for Kaggle docking workflows.

Creates standardized directory tree for receptors, ligands, outputs, and validation.
"""

from pathlib import Path
import subprocess


def create_docking_dirs(base_path="/kaggle/working/docking", subdirs=None):
    """Create standard docking directory structure.

    Args:
        base_path (str): Root directory for docking project. Default /kaggle/working/docking
        subdirs (list, optional): Subdirectory names. Defaults to standard structure.

    Returns:
        dict: {name: Path} for each created directory

    Example:
        >>> dirs = create_docking_dirs()
        >>> print(dirs['receptors'])  # /kaggle/working/docking/receptors
    """
    if subdirs is None:
        subdirs = ["receptors", "ligands", "outputs", "validation"]

    base = Path(base_path)
    created = {}

    for name in subdirs:
        dir_path = base / name
        dir_path.mkdir(parents=True, exist_ok=True)
        created[name] = dir_path
        print(f"  ✓ {dir_path}")

    return created


def show_directory_tree(base_path="/kaggle/working/docking"):
    """Display directory tree structure.

    Args:
        base_path (str): Root directory to display

    Example:
        >>> show_directory_tree()
        /kaggle/working/docking
        /kaggle/working/docking/ligands
        /kaggle/working/docking/outputs
        ...
    """
    print("\nLocal directory tree:")
    result = subprocess.run(
        ["find", base_path, "-type", "d"],
        capture_output=True, text=True, check=False
    )

    if result.returncode == 0:
        for line in sorted(result.stdout.strip().split('\n')):
            print(f"  {line}")
    else:
        print(f"  (unable to display tree: {result.stderr})")


def setup_kaggle_dirs(base_path="/kaggle/working/docking", show_tree=True):
    """Create directory structure and optionally display tree.

    Args:
        base_path (str): Root directory path
        show_tree (bool): Display directory tree after creation. Default True.

    Returns:
        dict: {name: Path} for each created directory

    Example:
        >>> dirs = setup_kaggle_dirs()
        >>> receptor_dir = dirs['receptors']
    """
    dirs = create_docking_dirs(base_path)

    if show_tree:
        show_directory_tree(base_path)

    return dirs
