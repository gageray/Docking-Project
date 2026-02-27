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
