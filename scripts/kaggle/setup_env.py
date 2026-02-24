"""Kaggle environment setup: install dependencies and verify tools.

Installs APT packages, Python packages, and GNINA binary.
"""

import subprocess
import importlib


def install_apt_packages(packages=None, quiet=True):
    """Install APT packages via apt-get.

    Args:
        packages (list, optional): Package names. Defaults to docking requirements.
        quiet (bool): Suppress verbose output. Default True.

    Example:
        >>> install_apt_packages(['openbabel', 'libopenbabel-dev'])
    """
    if packages is None:
        packages = ['openbabel', 'libopenbabel-dev']

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
        packages (list, optional): Package names. Defaults to docking requirements.
        break_system (bool): Use --break-system-packages flag. Default True.

    Example:
        >>> install_pip_packages(['rdkit', 'biopython'])
    """
    if packages is None:
        packages = [
            'pdbfixer',
            'openmm',
            'biopython',
            'mdanalysis',
            'openbabel-wheel',
            'rdkit'
        ]

    print("=== Installing Python packages ===")

    cmd = ["pip", "install"]
    if break_system:
        cmd.append("--break-system-packages")

    for pkg in packages:
        pkg_cmd = cmd + [pkg]
        subprocess.run(pkg_cmd, check=True, capture_output=True)
        print(f"  ✓ {pkg}")


def install_gnina(install_path="/usr/local/bin/gnina", force=False):
    """Download and install GNINA binary.

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

    print("=== Installing GNINA ===")

    if os.path.exists(install_path) and not force:
        print(f"  ✓ GNINA already present at {install_path}")
        return False

    # Download GNINA v1.1
    url = "https://github.com/gnina/gnina/releases/download/v1.1/gnina"
    subprocess.run(["wget", "--progress=bar", url, "-O", install_path], check=True)
    subprocess.run(["chmod", "+x", install_path], check=True)

    print(f"  ✓ GNINA installed to {install_path}")
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
    required_packages = [
        "rdkit", "openmm", "pdbfixer", "Bio", "openbabel",
        "MDAnalysis", "numpy", "pandas", "googleapiclient"
    ]

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

    if verify:
        return verify_environment()

    return None
