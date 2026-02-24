#!/bin/bash
# Usage:
#   bash kaggle_setup.sh           # compute-only (GNINA + analysis)
#   bash kaggle_setup.sh --full    # full stack (can re-run prep scripts too)

set -e

echo "=== Docking Repo Kaggle Setup ==="

# Always install the base compute dependencies
pip install -q -r requirements_kaggle.txt
echo "[OK] Base dependencies installed"

# If --full flag is passed, also install the heavy structural bio libs
if [[ "$1" == "--full" ]]; then
    echo "[INFO] Full build requested. Installing openmm, pdbfixer, pymol..."
    apt-get update -y -qq
    apt-get install -y -qq libopenmm-dev freeglut3-dev libxml2-dev 2>/dev/null
    pip install -q openmm pymol-open-source
    pip install -q git+https://github.com/openmm/pdbfixer.git
    echo "[OK] Full dependencies installed"
else
    echo "[INFO] Compute-only build. Skipping openmm/pdbfixer/pymol."
    echo "       Re-run with --full if you need to run receptor_prep or ligand_prep on Kaggle."
fi

echo "=== Setup Complete ==="
