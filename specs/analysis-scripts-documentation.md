# Analysis Scripts Documentation

**Purpose**: Post-docking analysis and validation tools for molecular docking pipeline.

**Location**: `scripts/analysis/`

**Data Output**: `data/analysis/`

---

## Overview

Three validation scripts for different stages of the docking pipeline:

1. **`compare_poses.py`** - Validate docked poses against crystal structure (no alignment)
2. **`compare_input_conformer.py`** - Validate input conformer against crystal structure (with alignment)
3. **`validate_conformer_generation.py`** - Validate conformer generation from SMILES
4. **`validate_conformer_generation_crest.py`** - Validate conformer generation using CREST/GFN2-xTB for charged ligands

---

## 1. compare_poses.py

### Purpose
Compare docked poses to reference crystal structure ligand in the same model space (NO alignment). Measures how well docking reproduced the crystal binding pose.

### Key Features
- Extracts reference ligand from receptor PDB using metadata
- Loads all docked poses from GNINA output SDF
- Calculates heavy-atom RMSD (hydrogens removed)
- Symmetry-aware RMSD calculation via RDKit
- Ranks poses by RMSD
- Generates JSON + CSV reports

### Usage
```bash
python scripts/analysis/compare_poses.py \
  --receptor data/receptors/6X3U_aligned.pdb \
  --metadata data/receptors/6x3u_metadata.json \
  --poses data/output/6X3U_aligned_apo_flumazenil_out.sdf \
  --output data/analysis/6X3U_pose_comparison
```

### Inputs
- **Receptor PDB**: Aligned receptor with reference ligand (HETATM records)
- **Metadata JSON**: Contains `target_ligand_resn` and `target_ligand_chain`
- **Poses SDF**: GNINA output with multiple docked poses
- **Output path**: Base path for reports (no extension)

### Outputs

#### JSON Report (`*_comparison.json`)
```json
{
  "receptor": "path/to/receptor.pdb",
  "reference_ligand": "FYP:Z",
  "docked_sdf": "path/to/poses.sdf",
  "total_poses": 20,
  "passing_poses": 5,
  "best_rmsd": 4.044,
  "results": [
    {
      "pose_idx": 0,
      "rmsd": 4.114,
      "gnina_score": -7.198,
      "cnn_score": 0.823,
      "cnn_affinity": -7.15,
      "pass_threshold": false
    }
  ]
}
```

#### CSV Report (`*_comparison.csv`)
```csv
pose_idx,rmsd,gnina_score,cnn_score,cnn_affinity,pass_threshold
0,4.1140,-7.1980,0.8230,-7.1500,False
1,4.1540,-7.0810,0.7950,-7.0800,False
```

#### Console Output
- Progress per pose (RMSD + GNINA score)
- Summary table with top 5 poses
- Pass/fail count (threshold: 2.0 Å)

### Implementation Notes
- **No alignment**: Poses and crystal ligand already in same coordinate system (aligned receptor)
- **Heavy atoms only**: `removeHs=True` in SDMolSupplier
- **Symmetry handling**: Uses `rdMolAlign.GetBestRMS()` for automorphism
- **GNINA properties**: Extracts `minimizedAffinity`, `CNNscore`, `CNNaffinity` from SDF

### Validation Threshold
- **RMSD < 2.0 Å**: Standard for successful pose reproduction
- **RMSD 2.0-3.0 Å**: Marginal
- **RMSD > 3.0 Å**: Failed to reproduce binding mode

---

## 2. compare_input_conformer.py

### Purpose
Compare input ligand conformer to reference crystal structure WITH alignment. Tests if the starting conformer geometry is similar to the crystal pose shape.

### Key Features
- Extracts reference ligand from receptor PDB
- Loads input ligand SDF (pre-docking)
- Symmetry-aware alignment and RMSD
- Tests conformational similarity (shape matching)
- Single conformer comparison

### Usage
```bash
python scripts/analysis/compare_input_conformer.py \
  --receptor data/receptors/6X3U_aligned.pdb \
  --metadata data/receptors/6x3u_metadata.json \
  --input data/ligands/flumazenil.sdf \
  --output data/analysis/6X3U_input_conformer_comparison
```

### Inputs
- **Receptor PDB**: Aligned receptor with reference ligand
- **Metadata JSON**: Ligand residue name and chain
- **Input SDF**: Pre-docking ligand conformer from ligand_prep.py
- **Output path**: Base path for reports

### Outputs

#### JSON Report (`*_comparison.json`)
```json
{
  "receptor": "path/to/receptor.pdb",
  "reference_ligand": "FYP:Z",
  "input_sdf": "path/to/input.sdf",
  "total_conformers": 1,
  "passing_conformers": 0,
  "best_rmsd": 2.470,
  "results": [
    {
      "conformer_idx": 0,
      "rmsd_aligned": 2.470,
      "pass_threshold": false
    }
  ]
}
```

#### CSV Report (`*_comparison.csv`)
```csv
conformer_idx,rmsd_aligned,pass_threshold
0,2.4700,False
```

#### Console Output
- Conformer count
- Aligned RMSD
- Pass/fail (threshold: 2.0 Å)
- Warning if starting conformer is far from crystal

### Implementation Notes
- **WITH alignment**: Uses `rdMolAlign.GetBestRMS()` which includes alignment
- **Shape comparison**: Tests if conformer geometry matches crystal independent of position
- **Heavy atoms only**: Protonation state doesn't affect RMSD
- **Single conformer expected**: Most ligand prep outputs 1 conformer

### Use Case
Diagnose whether poor docking results are due to:
- Bad starting conformer (high RMSD here)
- Docking failure (good RMSD here, bad in compare_poses.py)

---

## 3. validate_conformer_generation.py

### Purpose
Comprehensive validation of conformer generation from SMILES. Tests geometry quality, energy distribution, diversity, and crystal pose recovery.

### Key Features
- Generate N conformers from SMILES using ETKDG
- MMFF energy calculation and filtering
- Geometric validity checks (steric clashes)
- Conformational diversity analysis (pairwise RMSD)
- Crystal structure comparison (if provided)
- Energy distribution analysis
- Induced-fit binding detection

### Usage
```bash
python scripts/analysis/validate_conformer_generation.py \
  --smiles "CCOC(=O)c1ncn2c1C[N@@H+](C)C(=O)c1cc(F)ccc1-2" \
  --num-conformers 100 \
  --energy-window 40.0 \
  --receptor data/receptors/6X3U_aligned.pdb \
  --metadata data/receptors/6x3u_metadata.json \
  --output data/analysis/flumazenil_conformer_validation \
  --save-sdf data/ligands/flumazenil_validated_conformers.sdf
```

### Inputs
- **SMILES**: Input SMILES string (can be charged/protonated)
- **num-conformers**: Number of conformers to generate (default: 50)
- **energy-window**: Energy cutoff in kcal/mol above minimum (default: 40.0)
- **receptor**: Optional, for crystal comparison
- **metadata**: Optional, for crystal comparison
- **save-sdf**: Optional, output valid conformers to SDF

### Outputs

#### JSON Report (`*_validation.json`)
```json
{
  "smiles": "SMILES_STRING",
  "num_conformers_generated": 100,
  "num_conformers_valid": 1,
  "energy_window": 40.0,
  "min_energy": 2984542165.92,
  "geometry": {
    "num_atoms": 22,
    "num_bonds": 24,
    "num_rotatable_bonds": 2,
    "molecular_weight": 303.32,
    "conformers": [
      {
        "conf_id": 48,
        "energy": 2984542165.92,
        "min_nonbonded_dist": 2.15,
        "has_clash": false
      }
    ]
  },
  "diversity": {
    "num_conformers": 1,
    "mean_pairwise_rmsd": null,
    "median_pairwise_rmsd": null,
    "min_pairwise_rmsd": null,
    "max_pairwise_rmsd": null,
    "rmsd_matrix": []
  },
  "crystal_comparison": {
    "reference_ligand": "FYP:Z",
    "best_rmsd": 2.473,
    "best_conf_id": 48,
    "best_conf_energy_rank": 0,
    "best_conf_delta_energy": 0.0,
    "best_conf_energy_bin": "0-10",
    "passing_conformers": 0,
    "all_results": [
      {
        "conf_id": 48,
        "energy": 2984542165.92,
        "delta_energy": 0.0,
        "energy_bin": "0-10",
        "rmsd_to_crystal": 2.473,
        "pass_threshold": false
      }
    ]
  }
}
```

#### CSV Report (`*_validation.csv`)
```csv
conf_id,energy,delta_energy,energy_bin,min_nonbonded_dist,has_clash,rmsd_to_crystal,pass_threshold
48,2984542165.9200,0.0000,0-10,2.1500,False,2.4730,False
```

#### Console Output
```
Generating 100 conformers from SMILES: ...
Optimizing conformers with MMFF...

Energy Distribution:
  ΔE 0-10 kcal/mol: 1 conformers
  ΔE 100-inf kcal/mol: 99 conformers

Generated 100 conformers
Min energy: 2984542165.92 kcal/mol
Valid conformers (ΔE < 40.0 kcal/mol): 1

Checking geometric validity...
Conformers with steric clashes: 0/1

Checking conformational diversity...
WARNING: Only 1 conformer, cannot assess diversity

Comparing to crystal structure...
Best conformer: ID=48, RMSD=2.473 Å
Energy: ΔE=0.00 kcal/mol (bin: 0-10 kcal/mol)
Energy rank: 0/100
Conformers within 2.0 Å of crystal: 0/1
✗ Failed to recover crystal-like pose

============================================================
CONFORMER GENERATION VALIDATION SUMMARY
============================================================
SMILES: CCOC(=O)c1ncn2c1C[N@@H+](C)C(=O)c1cc(F)ccc1-2
Conformers generated: 100
Valid conformers: 1
Rotatable bonds: 2

Geometry:
  Steric clashes: 0/1

Diversity:
  Mean pairwise RMSD: N/A (only 1 conformer)

Crystal comparison:
  Best RMSD: 2.473 Å
  Energy rank: 0/100
  Passing conformers: 0/1
============================================================
```

### Validation Checks

#### 1. Geometric Validity
- **Steric clashes**: Atoms closer than 1.5 Å (non-bonded)
- **Min nonbonded distance**: Reports closest approach
- **Flags**: `has_clash` boolean per conformer

#### 2. Energy Distribution
- **Bins**: 0-10, 10-20, 20-30, 30-40, 40-50, 50-100, 100+ kcal/mol
- **Filtering**: Only conformers within `energy_window` are "valid"
- **Interpretation**: High bin counts at >40 kcal/mol indicate force field failures

#### 3. Conformational Diversity
- **Pairwise RMSD**: All-vs-all RMSD matrix
- **Statistics**: Mean, median, min, max pairwise RMSD
- **Threshold**: Mean RMSD < 1.0 Å indicates poor diversity
- **Note**: Requires ≥2 conformers

#### 4. Crystal Pose Recovery
- **Aligned RMSD**: Best RMS after alignment to crystal ligand
- **Energy rank**: Where the best-matching conformer ranks by energy
- **Energy bin**: Which energy window contains the best match
- **Induced-fit detection**: If best match is >20 kcal/mol above minimum

### Implementation Notes

#### Conformer Generation
```python
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.numThreads = 0
AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params)
```

#### Energy Calculation
```python
# MMFF optimization
AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0)

# Energy for each conformer
props = AllChem.MMFFGetMoleculeProperties(mol_h)
ff = AllChem.MMFFGetMoleculeForceField(mol_h, props, confId=conf_id)
energy = ff.CalcEnergy()
```

#### Energy Filtering
```python
min_energy = min(energies)
valid = [conf_id for conf_id, e in enumerate(energies)
         if e - min_energy <= energy_window]
```

#### Diversity Calculation
```python
for i, j in combinations:
    rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=i, refId=j)
    rmsd_matrix[i, j] = rmsd
```

### Energy Window Selection

| Window (kcal/mol) | Use Case |
|-------------------|----------|
| 10-20 | Standard thermal energy at 298K (~RT) |
| 20-30 | Moderately accessible conformations |
| 30-40 | Higher-energy states, potential induced-fit |
| 40+ | Non-physiological, likely force field failures |

**Default: 40.0 kcal/mol** - catches potential induced-fit conformations while filtering obvious failures.

### Induced-Fit Detection

If the best crystal-matching conformer has ΔE > 20 kcal/mol:
```
⚠ INDUCED-FIT BINDING DETECTED:
  Crystal pose is XX.X kcal/mol above solution minimum
  Binding site must provide ≥XX.X kcal/mol stabilization
```

Indicates the binding site stabilizes a higher-energy conformation through:
Indicates the binding site stabilizes a higher-energy conformation through:
- Hydrogen bonds (2-5 kcal/mol each)
- π-π stacking (2-4 kcal/mol)
- Electrostatic interactions (5-10 kcal/mol)
- Hydrophobic burial

---

## 4. validate_conformer_generation_crest.py

### Purpose
Validates conformer generation using CREST/GFN2-xTB, specifically designed for charged species where empirical force fields (like MMFF) fail.

### Key Features
- Generates a neutral 3D seed using ETKDG
- Performs conformational sampling with CREST at the GFN2-xTB level
- Uses fast search `--quick` modes to keep execution times practical
- Adds `--nocross` flag to disable genetic crossing for rigid molecules
- Computes energies (within a `--energy-window` default of 5.0) and compares the conformer ensemble against a crystal structure
- Extracts relevant geometries from CREST multi-XYZ outputs and maps them back to the charged SMILES topology

### Usage
```bash
python scripts/analysis/validate_conformer_generation_crest.py \
  --neutral-smiles "YOUR_NEUTRAL_SMILES" \
  --charged-smiles "YOUR_CHARGED_SMILES" \
  --energy-window 5.0 \
  --search-level quick \
  --receptor data/receptors/RECEPTOR.pdb \
  --metadata data/receptors/METADATA.json \
  --output data/analysis/LIGAND_crest_validation \
  --save-sdf data/ligands/LIGAND_crest_conformers.sdf \
  --threads 4
```

### Search Levels
- `default`: Full iMTD-GC search (extremely slow, thorough)
- `quick`: Reduced settings for a crude ensemble (recommended for high throughput)
- `squick`: Faster than `quick`
- `mquick`: Most aggressive reduction for maximum speed

---

## Workflow Integration

### Standard Validation Workflow

```bash
# 1. Validate conformer generation from SMILES
python scripts/analysis/validate_conformer_generation.py \
  --smiles "YOUR_SMILES" \
  --num-conformers 100 \
  --receptor data/receptors/RECEPTOR.pdb \
  --metadata data/receptors/METADATA.json \
  --output data/analysis/LIGAND_conformer_validation

# 2. Compare input conformer to crystal (shape similarity)
python scripts/analysis/compare_input_conformer.py \
  --receptor data/receptors/RECEPTOR.pdb \
  --metadata data/receptors/METADATA.json \
  --input data/ligands/LIGAND.sdf \
  --output data/analysis/LIGAND_input_comparison

# 3. Run docking (via GNINA)
# ... docking step ...

# 4. Compare docked poses to crystal (binding pose accuracy)
python scripts/analysis/compare_poses.py \
  --receptor data/receptors/RECEPTOR.pdb \
  --metadata data/receptors/METADATA.json \
  --poses data/output/DOCKED_POSES.sdf \
  --output data/analysis/LIGAND_pose_comparison
```

### Interpretation Matrix

| Step 1 Validation | Step 2 Input | Step 4 Poses | Diagnosis |
|-------------------|--------------|--------------|-----------|
| Many valid conformers | Low RMSD | Low RMSD | ✓ Pipeline working |
| Many valid conformers | Low RMSD | High RMSD | Docking failure (scoring/search) |
| Many valid conformers | High RMSD | High RMSD | Wrong conformer selected for docking |
| Few valid conformers | High RMSD | High RMSD | Conformer generation failure |
| Few valid conformers | Low RMSD | High RMSD | Docking failure (unlikely) |

---

## Dependencies

### Python Packages
- **RDKit**: Conformer generation, MMFF, RMSD
- **NumPy**: Energy/RMSD calculations
- **Biopython**: PDB parsing, ligand extraction
- **JSON/CSV**: Standard library for reports

### Conda Environment
```bash
conda activate docking
```

All scripts require the `docking` conda environment with RDKit, Biopython, NumPy installed.

---

## File Naming Conventions

### Analysis Outputs
- `{TARGET}_pose_comparison.json/csv` - Docked pose validation
- `{TARGET}_input_conformer_comparison.json/csv` - Input conformer validation
- `{LIGAND}_conformer_validation.json/csv` - Conformer generation validation
- `{LIGAND}_validated_conformers.sdf` - Valid conformers output

### Location
- **Scripts**: `scripts/analysis/`
- **Data**: `data/analysis/`
- **Validated conformers**: `data/ligands/`

---

## Common Issues

### Issue 1: Atom Count Mismatch
**Error**: `Atom count mismatch: generated=17, crystal=22`

**Cause**: Wrong SMILES or wrong crystal structure ligand

**Solution**: Verify SMILES matches crystal ligand, check metadata for correct `target_ligand_resn`

### Issue 2: 99% Energy Failures
**Output**: `ΔE 100-inf kcal/mol: 99 conformers`

**Cause**: Force field failure for charged/strained molecules

**Solution**:
- Test neutral SMILES
- Use alternative force field (PM6, GFN2-xTB)
- Use crystal structure geometry directly

### Issue 3: JSON Serialization Error
**Error**: `Object of type bool_ is not JSON serializable`

**Cause**: NumPy types not converted to Python native types

**Solution**: Already fixed - casts to `int()`, `float()`, `bool()` before JSON dump

### Issue 4: Single Conformer Diversity Check
**Output**: `WARNING: Only 1 conformer, cannot assess diversity`

**Cause**: Only 1 valid conformer found (expected for rigid molecules or force field failures)

**Solution**: Not an error - diversity check skipped, continues with validation

---

## Future Enhancements

### Potential Additions
1. **Interaction fingerprints**: Analyze H-bonds, hydrophobic contacts
2. **Torsion angle analysis**: Compare dihedral angles to crystal
3. **Binding affinity correlation**: Plot RMSD vs GNINA score
4. **Multi-receptor validation**: Batch validation across multiple targets
5. **PyMOL visualization scripts**: Auto-generate alignments for visual inspection

### Performance Optimizations
- Parallel conformer generation
- Cached crystal structure parsing
- Incremental validation (resume from checkpoint)

---

## References

### RDKit Documentation
- ETKDG: https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
- MMFF: https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html
- Alignment: https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html

### GNINA Output Format
- SDF properties: `minimizedAffinity`, `CNNscore`, `CNNaffinity`
- Multiple conformers per SDF file

### Validation Thresholds
- RMSD < 2.0 Å: Standard for pose reproduction (Ref: multiple docking benchmark papers)
- Energy window 20 kcal/mol: ~RT at 298K thermal energy
- Steric clash < 1.5 Å: VDW radii sum threshold
