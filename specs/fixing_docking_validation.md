# Fixing Docking Validation Issues

## Current Problem
**RMSD validation shows 2.5Å RMSD** - all poses are **POOR** quality and failed to reproduce the native crystal structure binding mode.

## Results Summary

### RMSD Results
```
Best Pose: #2 with 2.462 Å RMSD (POOR)
All 9 poses: 2.46-2.54 Å RMSD (all POOR)
```

### Docking Scores (ACTUAL PROBLEM IDENTIFIED)
```
Pose | minimizedAffinity | CNNscore | CNNaffinity | CNN_VS
-----|-------------------|----------|-------------|--------
  1  |       0.0         |  0.851   |   +4.04     | 3.44
  2  |       0.0         |  0.850   |   +3.97     | N/A
  3  |       0.0         |  0.850   |   +4.00     | N/A
  4  |       0.0         |  0.849   |   +4.03     | N/A
  ... (all poses similar)
```

**CRITICAL FINDINGS**:
1. **minimizedAffinity = 0.0** - Vina scores not saved (using `--cnn_scoring rescore` mode)
2. **CNNaffinity = +4 kcal/mol** - POSITIVE energy = BAD BINDING (should be negative)
3. **CNNscore = 0.85** - Decent confidence, but for a BAD pose
4. **All poses have similar bad scores** - Algorithm converged to wrong binding mode

**CONCLUSION**: Docking completely failed. GNINA thinks all poses are poor binders (positive energy) and RMSD confirms they're in the wrong location.

## Potential Issues & Fixes

### 1. **CNN Scoring Mode Issue** ⚠️ ROOT CAUSE
Current GNINA command from `scripts/kaggle/gnina_worker.py` line 58:
```python
"--cnn_scoring", "rescore"
```

**Problem**: `--cnn_scoring rescore` mode:
- Runs Vina docking to generate poses
- Then rescores with CNN
- **Does NOT save Vina energies** to output (that's why minimizedAffinity = 0.0)
- Only outputs CNN scores

**Actual scores we got**:
- `CNNaffinity: +4.04 kcal/mol` - **POSITIVE = BAD** (should be -8 to -11 for good binding)
- `CNNscore: 0.85` - Pose confidence (0-1), but confidence in a BAD pose
- `CNN_VS: 3.44` - Another CNN metric
- `CNNaffinity_variance: 0.53` - Uncertainty in prediction

**Fix**: Change scoring mode to get both scores:
```python
"--cnn_scoring", "all"  # Get both Vina AND CNN scores
# OR
"--cnn_scoring", "none"  # Just use Vina (faster, good baseline)
```

**Expected good scores**:
- Vina: -8 to -12 kcal/mol (NEGATIVE)
- CNNaffinity: -8 to -12 kcal/mol (NEGATIVE)
- CNNscore: > 0.9 for confident good poses

### 2. **Exhaustiveness Too Low** ⚠️ MAJOR ISSUE
Current run took ~60 seconds total - WAY too fast for accurate docking.

**Current setting**: `--exhaustiveness 8` (line 56 in gnina_worker.py)

**Problem**:
- Exhaustiveness=8 is bare minimum (quick screening)
- For validation/benchmarking, need exhaustiveness=16-32
- Higher exhaustiveness = more thorough search = better poses

**Fix - Update gnina_worker.py line 56**:
```python
"--exhaustiveness", "32",  # Increase from "8" to "32"
"--num_modes", "20",       # Increase from "9" to "20"
```

**Complete fixed command**:
```bash
gnina \
  -r receptor.pdbqt \
  -l ligand.sdf \
  -o output.sdf \
  --center_x X --center_y Y --center_z Z \
  --size_x SX --size_y SY --size_z SZ \
  --exhaustiveness 32 \
  --num_modes 20 \
  --cnn_scoring all \        # Changed from "rescore"
  --cpu 8
```

**Expected time**: 5-15 minutes per ligand with exhaustiveness=32 (vs 60 seconds with exhaustiveness=8)

### 3. **Input Conformer Mismatch**
We used `data/ligands/flumazenil.sdf` (RDKit-generated from SMILES) as reference.

**Problem**: This conformer might not match the crystal structure geometry
- Different ring puckers
- Different torsion angles
- Different protonation state

**Fix - Compare conformers**:
```python
# Compare our generated conformer vs crystal structure
python scripts/pipeline/score_native_rmsd.py \
  --native data/ligands/flumazenil.sdf \
  --docked data/output/native_FYP_crystal_CORRECT.sdf \
  --out conformer_comparison.txt
```

If RMSD > 0.5Å, our input conformer is wrong.

**Fix - Use crystal conformer as input**:
Instead of generating from SMILES, extract FYP from crystal and use that as docking input.

### 4. **Box Definition Issues**
Binding box might not cover the full binding site.

**Check**:
- View box in PyMOL with `box_params.json`
- Ensure box includes His102, Phe77, and all key residues
- Box should be centered on native ligand position

**Fix**: Adjust box size/center in `scripts/receptor/box_definition.py`

### 5. **Receptor Preparation Issues**
- Missing hydrogens
- Wrong protonation states
- Incorrect side chain conformations

**Check**:
```bash
# Verify receptor has hydrogens
grep "^ATOM.*H" data/receptors/prepped/6X3U_apo.pdbqt | wc -l
```

### 6. **Ligand Protonation State**
Flumazenil might need different protonation at pH 7.4

**Fix**:
- Check if Dimorphite-DL was used
- Verify nitrogen protonation states match physiological pH

## Step-by-Step Debugging Protocol

### Step 1: Check Docking Scores
```bash
conda activate docking
python scripts/pipeline/extract_scores.py
# Look for:
# - minimizedAffinity < -8 kcal/mol (good binding)
# - CNNscore > 0.8 (high confidence)
```

### Step 2: Check Conformer Match
```python
# Does our input match the crystal?
from rdkit.Chem import rdMolAlign
rmsd_input_vs_crystal = compare_conformers()
# If > 0.5Å, input conformer is wrong
```

### Step 3: Re-run with Higher Exhaustiveness
```bash
# Use proper GNINA parameters
gnina --receptor ... --ligand ... \
  --exhaustiveness 32 \
  --num_modes 20 \
  --cpu 8
```

### Step 4: Visual Inspection
```bash
pymol data/receptors/6X3U_aligned.pdb \
     data/output/6X3U_apo_flumazenil_out.sdf
# Manually check if ANY pose looks correct
```

## Expected Benchmarks

### Good Self-Docking Performance:
- **RMSD < 2.0 Å**: Acceptable
- **RMSD < 1.0 Å**: Excellent
- **Best pose in top 3**: Good ranking function
- **minimizedAffinity**: -9 to -11 kcal/mol (typical benzodiazepine)

### Our Current Performance:
- **RMSD 2.5 Å**: FAILED ❌
- **All poses similar**: Converged to wrong binding mode ❌
- **Fast runtime (60s)**: Insufficient search ❌

## Files to Check
- Docking output: `data/output/6X3U_apo_flumazenil_out.sdf`
- Scores: Extract with `extract_scores.py`
- Native reference: `data/output/native_FYP_crystal_CORRECT.sdf`
- RMSD results: `data/output/flumazenil_rmsd_CORRECT.txt`

## Next Steps - Priority Order

### IMMEDIATE FIXES (Do These First)
1. **Update gnina_worker.py**:
   - Line 56: Change `"--exhaustiveness", "8"` to `"--exhaustiveness", "32"`
   - Line 57: Change `"--num_modes", "9"` to `"--num_modes", "20"`
   - Line 58: Change `"--cnn_scoring", "rescore"` to `"--cnn_scoring", "all"`

2. **Re-run flumazenil docking** with fixed parameters

3. **Extract scores** and verify:
   ```bash
   # Check that we now get NEGATIVE energies
   python -c "from rdkit import Chem; mol = Chem.SDMolSupplier('output.sdf')[0]; print(mol.GetProp('minimizedAffinity'))"
   # Should be -8 to -11 kcal/mol (NEGATIVE)
   ```

4. **Re-run RMSD validation** with new docking results

### SECONDARY CHECKS (If Still Failing)
5. **Compare input conformer to crystal** structure conformer
6. **Verify box definition** covers binding site properly
7. **Check receptor prep** for hydrogens and protonation

## Success Criteria
- At least 1 pose with RMSD < 2.0 Å
- That pose should be in top 5 by score
- Exhaustiveness run should take 5-15 minutes
- Visual inspection shows correct binding mode

---

## Summary

**ROOT CAUSES IDENTIFIED**:
1. ⚠️ **Exhaustiveness too low (8)** - Insufficient conformational search
2. ⚠️ **CNN scoring mode wrong** - Using "rescore" instead of "all", losing Vina scores
3. ⚠️ **Positive CNNaffinity (+4 kcal/mol)** - GNINA itself thinks these are bad poses

**EVIDENCE**:
- All poses: 2.5Å RMSD (POOR) ❌
- All poses: +4 kcal/mol CNNaffinity (BAD BINDING) ❌
- Runtime: 60 seconds (TOO FAST) ❌
- minimizedAffinity: 0.0 (NOT SAVED) ❌

**STATUS**: Docking validation COMPLETELY FAILED. Must fix parameters and re-run before using for any novel compounds.

---

## RESOLUTION - IMPLEMENTATION COMPLETE

### Code Changes Implemented

#### 1. Created `scripts/pipeline/gnina_runner.py`
Centralized GNINA execution logic with configuration profiles:
- **GNINAConfig class**: Immutable configuration with 4 profiles
  - `quick_screening`: exhaustiveness=8, num_modes=9, cnn_scoring=none
  - `standard`: exhaustiveness=16, num_modes=20, cnn_scoring=all
  - `thorough`: exhaustiveness=32, num_modes=20, cnn_scoring=all
  - `validation`: exhaustiveness=32, num_modes=20, cnn_scoring=all
- **GNINARunner class**: Execution logic with proper error handling
- **Score extraction**: Built-in parsing of all GNINA output scores
- **Validation**: Automatic quality checks on output files

#### 2. Refactored `scripts/kaggle/gnina_worker.py` to `scripts/kaggle/worker.py`
Generic multi-GPU task worker (GNINA-agnostic):
- Imports GNINARunner from pipeline module
- Loads configuration from work_queue.json metadata
- Supports multiple task types (extensible architecture)
- Proper score logging and validation warnings
- Backup saved as `gnina_worker_deprecated.py`

#### 3. Updated `scripts/pipeline/local_dispatcher.py`
Added metadata support for new worker format:
- New parameters: `gnina_profile`, `task_type`, `max_workers`
- Returns dict with metadata + jobs array
- CLI arguments for profile selection
- Updated job queue format

#### 4. Enhanced `scripts/pipeline/extract_scores.py`
Quality assessment and warning flags:
- `classify_binding()`: EXCELLENT/GOOD/MODERATE/WEAK/POOR_POSITIVE_ENERGY
- `check_warnings()`: VINA_NOT_SAVED, POSITIVE_CNN_ENERGY, LOW_CNN_CONFIDENCE
- New CSV columns: Binding_Quality, Warning

#### 5. Created `scripts/pipeline/run_validation.py`
Automated validation workflow:
- Runs GNINA with validation-grade parameters
- Calculates RMSD vs native pose
- Extracts and validates scores
- Generates markdown report with pass/fail criteria
- Success criteria:
  - RMSD < 2.0A
  - Best RMSD pose in top 5 by score
  - At least one pose with negative energy
  - Runtime > 300s (sufficient exhaustiveness)

#### 6. Updated `config.json`
Added validation and production configuration sections:
- `docking_validation`: Test cases, GNINA profile, output directory
- `docking_production`: Production settings, thorough profile

#### 7. Updated `scripts/pipeline/local_dispatcher.json`
Added new configuration parameters:
- `gnina_profile`: "thorough"
- `task_type`: "gnina_docking"
- `max_workers`: 2

#### 8. Created `specs/gnina_configuration.md`
Comprehensive GNINA parameter documentation:
- Profile comparison table
- Parameter explanations
- Score interpretation guide
- Common issues and fixes
- Usage examples
- Migration guide

### Configuration Changes (CRITICAL FIXES)

**OLD (BROKEN)**:
```python
"--exhaustiveness", "8"        # Too low
"--num_modes", "9"             # Too few
"--cnn_scoring", "rescore"     # Loses Vina scores
```

**NEW (FIXED)**:
```python
exhaustiveness: 32             # Thorough search
num_modes: 20                  # More diverse poses
cnn_scoring: "all"             # Both Vina AND CNN scores
```

### Validation Workflow

**Run before production**:
```bash
python scripts/pipeline/run_validation.py \
    --receptor 6X3U_apo \
    --ligand flumazenil \
    --native data/ligands/flumazenil.sdf \
    --config config.json \
    --out validation_report.md
```

**Expected Results**:
- RMSD < 2.0A (PASS)
- Negative energies in output (PASS)
- Runtime > 5 minutes (PASS)
- Best pose in top 5 (PASS)

### File Structure Summary

**NEW FILES**:
- `scripts/pipeline/gnina_runner.py` - GNINA execution module
- `scripts/pipeline/run_validation.py` - Automated validation
- `specs/gnina_configuration.md` - Documentation

**RENAMED FILES**:
- `scripts/kaggle/gnina_worker.py` → `scripts/kaggle/worker.py`
- Backup: `scripts/kaggle/gnina_worker_deprecated.py`

**MODIFIED FILES**:
- `scripts/pipeline/local_dispatcher.py` - Metadata support
- `scripts/pipeline/extract_scores.py` - Quality assessment
- `scripts/pipeline/local_dispatcher.json` - New parameters
- `config.json` - Validation/production sections
- `specs/fixing_docking_validation.md` - This resolution section

### Testing Before Production

1. **Run validation workflow**:
   ```bash
   cd /path/to/project
   python scripts/pipeline/run_validation.py \
       --receptor 6X3U_apo \
       --ligand flumazenil \
       --native data/ligands/flumazenil.sdf \
       --config config.json
   ```

2. **Review validation report**: Check for PASS status and RMSD < 2.0A

3. **Dispatch to Kaggle**:
   ```bash
   python scripts/pipeline/local_dispatcher.py \
       --config scripts/pipeline/local_dispatcher.json \
       --gnina-profile thorough
   ```

4. **Extract and analyze results**:
   ```bash
   python scripts/pipeline/extract_scores.py
   # Review Binding_Quality and Warning columns in CSV
   ```

### Architecture Benefits

**Separation of Concerns**:
- GNINA logic in pipeline module (reusable)
- Worker is task-agnostic (extensible)
- Configuration-driven behavior (maintainable)

**Quality Assurance**:
- Automated validation before production
- Quality flags in score extraction
- Proper error handling and logging

**Maintainability**:
- Single source of truth for GNINA parameters
- No hardcoded values in worker
- Comprehensive documentation

### Migration Notes

**For Kaggle Kernels**:
1. Update kernel to use `worker.py` instead of `gnina_worker.py`
2. Ensure `scripts/pipeline/gnina_runner.py` is uploaded with kernel
3. Update `work_queue.json` to include metadata section

**For Local Testing**:
1. Use `gnina_runner.py` directly for single docking runs
2. Use `run_validation.py` for benchmarking
3. Use `local_dispatcher.py` for batch job creation

### Known Limitations

1. **Runtime**: Thorough profile (exhaustiveness=32) takes 10-20 min/ligand
   - This is NORMAL and expected for validation-grade docking
   - Kaggle kernels can handle 30-40 ligands per run (9-12 hour limit)

2. **Deprecated file kept**: `gnina_worker_deprecated.py` retained as backup
   - Can be deleted after confirming new worker functions correctly

3. **Box parameters**: Still loaded from `box_params.json`
   - Future: Integrate box calculation into validation workflow

### Next Steps

1. **Test validation locally** with flumazenil
2. **Review validation report** for PASS/FAIL
3. **Push updated worker to Kaggle** if validation passes
4. **Run production docking** with thorough profile
5. **Analyze results** with enhanced extract_scores.py

**Status**: All code changes implemented. Ready for validation testing.
