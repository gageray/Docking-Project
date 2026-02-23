# Molecular Docking Pipeline

## Overview
Automated molecular docking pipeline for GNINA. Processes receptors (CIF→PDB) and ligands (SMILES→SDF) with standardized geometry alignment.

## Stack
- **RDKit** - Ligand conformer generation
- **Dimorphite-DL** - pH 7.4 protonation
- **Biopython** - Structure parsing/alignment
- **OpenMM/PDBFixer** - Protein prep
- **PyMOL** - Visualization
- **NumPy** - Geometry transforms

## Pipeline
1. **Ligand Prep** (`scripts/ligand/`) - SMILES→SDF, conformer generation, grid search
2. **Receptor Prep** (`scripts/receptor/`) - CIF→PDB, chain extraction, H addition
3. **Geometry Align** (`scripts/receptor/align_model_space.py`) - Standardize to pocket-centered coords
4. **Box Definition** (`scripts/receptor/box_definition.py`) - Calculate binding boxes
5. **Visualization** (`scripts/visualization/`) - PyMOL validation/rendering

## Config
`config.json` = single source of truth. No hardcoded paths/params.

## Rules
- **NO HALTS** - Log errors, continue execution
- Use `utils.py` for shared logging
- pH 7.4, MMFF force field
- Energy filter: >20 kcal/mol above min
- Sequence-aware CA alignment

## Memory Bank Protocol
**ALWAYS read at session start:**
- `memory-bank/active-context.md` - Current state
- `memory-bank/progress.md` - Completed tasks
- `memory-bank/project-brief.md` - Overview

**Update at session end** with new work.

## Structure
```
scripts/
  ligand/          - Ligand prep
  receptor/        - Receptor prep
  visualization/   - PyMOL tools
specs/             - Design docs
data/              - Raw/processed (gitignored)
config.json        - Pipeline config
```
