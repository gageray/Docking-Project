# Project Brief

**Molecular Docking Pipeline for GNINA**

Automated data pipeline for preparing receptor-ligand complexes for molecular docking studies using GNINA software.

## Core Objectives
1. Convert protein structures (CIF → PDB) with chain extraction and repair
2. Generate ligand conformers (SMILES → SDF) with pH-dependent protonation
3. Standardize receptor geometry to pocket-centered coordinate system
4. Calculate binding box parameters via sequence-aware alignment
5. Provide PyMOL-based visualization for validation

## Key Features
- **Config-driven**: Single source of truth in `config.json`
- **Robust**: No-halt execution with comprehensive logging
- **Scientific accuracy**: pH 7.4, MMFF/Amber force fields, energy filtering
- **Modular**: Separated ligand/receptor/visualization scripts
- **Validated**: PyMOL scripts for visual checkpoint verification

## Target Receptors
- GABA_A receptor subtypes (alpha1, alpha2)
- Pentameric ligand-gated ion channels

## Target Ligands
- Diazepam, Flumazenil (single conformers)
- Ferulic acid (single conformer)
- Docosanyl ferulate (72-conformer grid search with locked tail)
