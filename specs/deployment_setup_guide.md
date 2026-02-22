# Docking Pipeline — Environment Setup Guide

## Architecture Overview

Two environments, one repo. You develop and validate locally in WSL, then clone the same repo to Kaggle for GPU-accelerated docking (GNINA on T4).

```
┌─────────────────────────────────┐     ┌──────────────────────────────┐
│  LOCAL (WSL2 Ubuntu)            │     │  REMOTE (Kaggle T4 GPU)      │
│                                 │     │                              │
│  Full stack:                    │     │  Compute-only:               │
│  - ligand_prep.py               │     │  - GNINA docking             │
│  - receptor_prep.py             │     │  - (optional) receptor_prep  │
│  - box_definition.py            │     │  - (optional) ligand_prep    │
│  - PyMOL visualization          │     │                              │
│  - All validation/debugging     │     │  No PyMOL, no vis scripts    │
│                                 │     │                              │
│  conda env (environment.yml)    │     │  pip install (kaggle_setup)  │
└──────────┬──────────────────────┘     └──────────┬───────────────────┘
           │                                       │
           └───────── same git repo ───────────────┘
```

---

## Part 1: Repository Structure

```
docking-repo/
├── config.json
├── environment.yml          # Local WSL conda env
├── requirements_kaggle.txt  # Kaggle pip dependencies
├── kaggle_setup.sh          # Kaggle bootstrap (with --full flag)
├── README.md
│
├── scripts/
│   ├── ligand_prep.py
│   ├── receptor_prep.py
│   ├── box_definition.py
│   ├── vis_stage1_ligands.py
│   ├── vis_stage2_box.py
│   └── vis_stage3_results.py
│
└── data/
    ├── receptors/
    │   ├── 6X3U.cif
    │   ├── 9CRS.cif
    │   └── 9CSB.cif
    └── ligands/
        └── (generated SDF files go here)
```

---

## Part 2: Local WSL Setup

### 2.1 — Install Miniforge (if you don't have conda yet)

```bash
# Inside WSL terminal
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh
# Accept defaults. Restart terminal after.
```

### 2.2 — environment.yml

```yaml
name: docking
channels:
  - conda-forge
dependencies:
  - python=3.10
  - openmm
  - pdbfixer
  - pymol-open-source
  - rdkit
  - biopython
  - numpy
  - pip
  - pip:
    - dimorphite-dl
```

### 2.3 — Create and activate

```bash
cd ~/docking-repo        # or wherever your repo lives
conda env create -f environment.yml
conda activate docking
```

### 2.4 — Verify installation

```bash
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "from pdbfixer import PDBFixer; print('PDBFixer OK')"
python -c "from openmm.app import PDBFile; print('OpenMM OK')"
python -c "import pymol; print('PyMOL OK')"
python -c "from Bio.PDB import MMCIFParser; print('BioPython OK')"
python -c "from dimorphite_dl import DimorphiteDL; print('Dimorphite OK')"
```

If any of these fail, `conda install -c conda-forge <package>` individually.

### 2.5 — Run the pipeline locally

```bash
conda activate docking
cd ~/docking-repo

python scripts/ligand_prep.py      # generates SDF files in data/ligands/
python scripts/receptor_prep.py    # generates prepped PDBs in data/receptors/
python scripts/box_definition.py   # generates box_params.json
```

### 2.6 — Visualize locally (PyMOL)

```bash
pymol
# Inside PyMOL console:
# run scripts/vis_stage1_ligands.py
# run scripts/vis_stage2_box.py
```

---

## Part 3: Kaggle Deployment

### 3.1 — requirements_kaggle.txt

Separate from the local environment.yml. These are the pip-installable packages for Kaggle's Docker base.

```
rdkit-pypi
dimorphite-dl
biopython
numpy
```

Note: `openmm`, `pdbfixer`, and `pymol-open-source` are intentionally excluded from the default Kaggle requirements. They're only needed if you want to re-run receptor_prep or ligand_prep on Kaggle (see the `--full` flag below).

### 3.2 — kaggle_setup.sh

```bash
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
```

### 3.3 — Kaggle Notebook Cells

```python
# Cell 1: Clone and setup
!git clone https://github.com/<you>/docking-repo.git /kaggle/working/docking-repo
%cd /kaggle/working/docking-repo

# Compute-only (most common — you already have prepped files from local)
!bash kaggle_setup.sh

# OR full build (if you want to re-run everything from scratch on Kaggle)
# !bash kaggle_setup.sh --full
```

```python
# Cell 2: Upload your locally-prepped data (if not committed to repo)
# Option A: commit prepped PDBs and SDFs to the repo (easiest)
# Option B: upload as a Kaggle Dataset and symlink
# !ln -s /kaggle/input/your-dataset/data ./data
```

```python
# Cell 3: Run GNINA docking
# (Your GNINA execution commands here)
```

### 3.4 — What to commit vs. what to upload as Kaggle Dataset

**Commit to repo** (small, text-like, changes rarely):
- All Python scripts
- config.json, environment.yml, requirements_kaggle.txt, kaggle_setup.sh
- Raw CIF files (a few MB each — fine for git)

**Upload as Kaggle Dataset** (large, binary, or generated):
- Prepped PDB files (if they're large after hydrogen addition)
- Generated SDF conformer libraries
- GNINA binary itself (if not already available on Kaggle)
- box_params.json (small, but it's a generated artifact)

Alternatively, if your total repo stays under ~50MB including data, just commit everything. Simpler.

---

## Part 4: The Editor / IDE Situation

### The Problem
Antigravity (your VS Code fork) doesn't support the WSL extension, which normally lets you edit WSL files with native Linux intellisense and execution.

### Option A: Work entirely in the WSL terminal (Recommended for this project)

This is what most people actually do for computational science pipelines. You don't need an IDE for this project — the scripts are straightforward, you're not building a web app.

```bash
# Open WSL terminal (Windows Terminal → Ubuntu tab)
cd ~/docking-repo
conda activate docking

# Edit with nano, vim, or micro (a modern terminal editor)
# Install micro if you want something user-friendly:
curl https://getmic.ro | bash
sudo mv micro /usr/local/bin/

# Then just:
micro scripts/ligand_prep.py
```

Run scripts, check output, edit, repeat. PyMOL gives you the visual feedback loop.

### Option B: SSH from Antigravity using Remote-SSH extension

If you really want GUI editing with syntax highlighting:

```bash
# Inside WSL, install and start SSH server
sudo apt install openssh-server
sudo service ssh start

# Check it's running
sudo service ssh status
```

Then in Antigravity, install the **Remote - SSH** extension (this is separate from the WSL extension and should work in any VS Code fork). Connect to:

```
ssh username@localhost
```

or

```
ssh username@$(hostname -I | awk '{print $1}')
```

This gives you the full remote editing experience — file explorer, integrated terminal, Python extension all running inside WSL. Functionally identical to the WSL extension.

**To auto-start SSH on WSL boot** (so you don't have to run it every time):

```bash
# Add to ~/.bashrc or create /etc/wsl.conf
echo '[boot]' | sudo tee -a /etc/wsl.conf
echo 'command=service ssh start' | sudo tee -a /etc/wsl.conf
```

### Option C: Install Antigravity natively inside WSL

This is possible but weird. You'd be running a GUI app from WSL via WSLg (Windows 11's built-in X server). It works, but:

- It'll feel slightly sluggish (GUI rendering over WSLg)
- You lose the Windows-native window management
- Extension ecosystem may have Linux compatibility issues
- You're now maintaining two Antigravity installs

**Verdict:** Don't bother. Option A or B are cleaner. If you need an IDE, SSH into WSL from the Windows-side Antigravity. If you don't (and for this project you probably don't), just use the terminal.

---

## Part 5: Git Workflow

### Initial setup (local WSL)

```bash
cd ~/docking-repo
git init
git remote add origin https://github.com/<you>/docking-repo.git

# .gitignore
cat > .gitignore << 'EOF'
__pycache__/
*.pyc
.agent/
memory-bank/
*.egg-info/

# Large generated files (if you'd rather upload as Kaggle Dataset)
# data/receptors/*_prepped.pdb
# data/ligands/*.sdf
# box_params.json
EOF

git add .
git commit -m "initial pipeline code"
git push -u origin main
```

### Iterate cycle

```
Local WSL                           Kaggle
─────────                           ──────
1. Edit scripts
2. Run pipeline, validate in PyMOL
3. git commit + push
                                    4. !git pull (or re-clone)
                                    5. !bash kaggle_setup.sh
                                    6. Run GNINA
                                    7. Download results
8. Load results in PyMOL locally
   (vis_stage3_results.py)
```

---

## Part 6: Checklist Before First Kaggle Run

- [ ] All scripts run locally without errors
- [ ] `box_params.json` generated and visually validated in PyMOL (Stage 2)
- [ ] Conformer SDFs validated in PyMOL (Stage 1)
- [ ] Prepped PDBs are in the repo or uploaded as Kaggle Dataset
- [ ] `kaggle_setup.sh` tested with both default and `--full` flags
- [ ] GNINA binary available on Kaggle (check if pre-installed or upload)
- [ ] `.gitignore` is clean — no massive binaries accidentally committed
