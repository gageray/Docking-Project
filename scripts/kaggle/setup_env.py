"""Kaggle environment setup: install dependencies and verify tools.

Installs APT packages, Python packages, and GNINA binary.
"""

import subprocess
import importlib


def install_apt_packages(packages=None, quiet=True):
    """Install APT packages via apt-get.

    Args:
        packages (list, optional): Package names. Defaults to None (skip APT install).
        quiet (bool): Suppress verbose output. Default True.

    Example:
        >>> install_apt_packages(['openbabel', 'libopenbabel-dev'])
    """
    if packages is None:
        print("=== Skipping APT packages (not needed for GNINA) ===")
        return

    print("=== Installing APT packages ===")

    # Update package list
    update_cmd = ["apt-get", "update", "-qq"] if quiet else ["apt-get", "update"]
    subprocess.run(update_cmd, check=False)

    # Install packages
    install_cmd = ["apt-get", "install", "-y"]
    if quiet:
        install_cmd.append("-qq")
    install_cmd.extend(packages)

    subprocess.run(install_cmd, check=True)
    print(f"  ✓ Installed {len(packages)} APT packages")


def install_pip_packages(packages=None, break_system=True):
    """Install Python packages via pip.

    Args:
        packages (list, optional): Package names. Defaults to RDKit only.
        break_system (bool): Use --break-system-packages flag. Default True.

    Example:
        >>> install_pip_packages(['rdkit'])
    """
    if packages is None:
        packages = ['rdkit']

    print("=== Installing Python packages ===")

    cmd = ["pip", "install"]
    if break_system:
        cmd.append("--break-system-packages")

    for pkg in packages:
        pkg_cmd = cmd + [pkg]
        subprocess.run(pkg_cmd, check=True, capture_output=True)
        print(f"  ✓ {pkg}")


def install_gnina(install_path="/usr/local/bin/gnina", force=False):
    """Install GNINA from Kaggle dataset.

    Args:
        install_path (str): Where to install GNINA. Default /usr/local/bin/gnina
        force (bool): Reinstall even if already present. Default False.

    Returns:
        bool: True if newly installed, False if already present

    Example:
        >>> install_gnina()
        >>> # GNINA now available as 'gnina' command
    """
    import os
    import shutil

    print("=== Installing GNINA ===")

    dataset_path = "/kaggle/input/gnina-bin/GNINA"

    if not os.path.exists(dataset_path):
        raise FileNotFoundError(f"GNINA dataset not found at {dataset_path}. Ensure 'ineptrobot/gnina-bin' is in dataset_sources.")

    print(f"  ✓ Found GNINA dataset at {dataset_path}")

    shutil.copy(dataset_path, install_path)
    subprocess.run(["chmod", "+x", install_path], check=True)
    print(f"  ✓ GNINA installed from dataset to {install_path}")
    return True


def verify_environment():
    """Verify GPU, GNINA, and Python packages are available.

    Returns:
        dict: Status of GPU, GNINA, and each Python package

    Example:
        >>> status = verify_environment()
        >>> if not status['packages']['rdkit']['available']:
        ...     print("RDKit missing!")
    """
    print("\n=== Verification ===")

    status = {
        "gpu": {},
        "gnina": {},
        "packages": {}
    }

    # Check GPU
    print("\n--- GPU ---")
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True, text=True, check=True
        )
        gpu_info = result.stdout.strip()
        print(f"  {gpu_info}")
        status["gpu"]["available"] = True
        status["gpu"]["info"] = gpu_info
    except Exception as e:
        print(f"  ✗ GPU check failed: {e}")
        status["gpu"]["available"] = False
        status["gpu"]["error"] = str(e)

    # Check GNINA
    print("\n--- GNINA ---")
    try:
        result = subprocess.run(
            ["gnina", "--version"],
            capture_output=True, text=True, check=False
        )
        version_line = result.stderr.split('\n')[0] if result.stderr else "installed"
        print(f"  ✓ {version_line}")
        status["gnina"]["available"] = True
        status["gnina"]["version"] = version_line
    except Exception as e:
        print(f"  ✗ GNINA not found: {e}")
        status["gnina"]["available"] = False
        status["gnina"]["error"] = str(e)

    # Check Python packages
    print("\n--- Python packages ---")
    required_packages = ["rdkit", "googleapiclient"]

    for pkg in required_packages:
        try:
            m = importlib.import_module(pkg)
            ver = getattr(m, "__version__", "ok")
            print(f"  ✓ {pkg:20s} {ver}")
            status["packages"][pkg] = {"available": True, "version": ver}
        except ImportError as e:
            print(f"  ✗ {pkg:20s} MISSING — {e}")
            status["packages"][pkg] = {"available": False, "error": str(e)}

    return status


def extract_ligand_zips(data_dir="/kaggle/working/ligands"):
    """Extract PDBQT files from zip archives in the ligand directory."""
    import zipfile
    import os
    from pathlib import Path

    p = Path(data_dir)
    if not p.exists():
        return

    print("\n=== Extracting Ligand Ensembles ===")
    for z_path in p.glob("*.zip"):
        try:
            with zipfile.ZipFile(z_path, 'r') as zf:
                # Only extract PDBQT files
                pdbqt_files = [name for name in zf.namelist() if name.endswith('.pdbqt')]
                for pdbqt_name in pdbqt_files:
                    zf.extract(pdbqt_name, p)
                print(f"  ✓ Extracted {len(pdbqt_files)} PDBQT files from {z_path.name}")
        except Exception as e:
            print(f"  ✗ Failed to extract {z_path.name}: {e}")


def setup_kaggle_environment(verify=True):
    """Complete Kaggle environment setup: install everything and verify.

    Args:
        verify (bool): Run verification after installation. Default True.

    Returns:
        dict or None: Verification status if verify=True, else None

    Example:
        >>> status = setup_kaggle_environment()
        >>> if status['gnina']['available']:
        ...     print("Ready to dock!")
    """
    install_apt_packages()
    install_pip_packages()
    install_gnina()
    extract_ligand_zips()

    if verify:
        return verify_environment()

    return None
