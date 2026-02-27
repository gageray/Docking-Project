# GNINA Configuration Profiles

## Profile Comparison

| Profile | Exhaustiveness | Num Modes | CNN Scoring | Use Case | Time/Ligand |
|---------|----------------|-----------|-------------|----------|-------------|
| quick_screening | 8 | 9 | none | Initial filter | 1-2 min |
| standard | 16 | 20 | all | Production runs | 5-10 min |
| thorough | 32 | 20 | all | Final candidates | 10-20 min |
| validation | 32 | 20 | all | Benchmarking | 10-20 min |

## Parameter Details

### exhaustiveness
Controls thoroughness of conformational search.
- **8**: Bare minimum (quick screening, NOT for validation)
- **16**: Standard production (good balance)
- **32**: Thorough search (validation-grade)
- **64**: Extremely thorough (research publications)

**Rule of thumb**: Higher = more accurate, but slower (roughly linear scaling)

### num_modes
Number of top-scoring poses to return.
- **9**: Quick screening
- **20**: Standard (captures diversity)
- **50**: Exhaustive analysis

### cnn_scoring
Which scoring functions to use:
- **none**: Vina only (fastest, good baseline)
- **rescore**: Vina search, then CNN rescore (LOSES Vina scores in output)
- **all**: Both Vina AND CNN (RECOMMENDED)

**CRITICAL**: Use `all` to get both minimizedAffinity AND CNNaffinity in output.

## Score Interpretation

### minimizedAffinity (Vina)
- < -10 kcal/mol: Excellent binding
- -8 to -10 kcal/mol: Good binding
- -6 to -8 kcal/mol: Moderate binding
- > -6 kcal/mol: Weak binding
- **> 0 kcal/mol**: BAD - indicates docking failure

### CNNaffinity
- Same scale as Vina
- Should correlate with Vina
- More accurate for novel scaffolds

### CNNscore
- Confidence: 0 (low) to 1 (high)
- > 0.9: High confidence
- 0.7-0.9: Moderate confidence
- < 0.7: Low confidence

## Validation Protocol

For self-docking benchmarks:
1. Use `validation` profile (exhaustiveness=32)
2. Calculate RMSD vs native pose
3. Check RMSD < 2.0A
4. Verify negative binding energies
5. Runtime should be 5-15 minutes (indicates sufficient search)

## Common Issues

### All Energies Positive
- Indicates docking failure
- Check box definition covers binding site
- Verify receptor preparation (protonation state)
- Try higher exhaustiveness

### Missing Vina Scores
- Using `cnn_scoring="rescore"` instead of "all"
- Change to `cnn_scoring="all"` to preserve both score types

### RMSD > 2.0A
- Insufficient exhaustiveness
- Box too small or misplaced
- Receptor preparation issues
- Try validation profile (exhaustiveness=32)

### Runtime < 5 minutes
- Exhaustiveness too low
- Not thorough enough for validation
- Increase to 32 for proper benchmarking

## Usage Examples

### Local Testing
```bash
python scripts/pipeline/gnina_runner.py \
    --receptor data/receptors/prepped/6X3U_apo.pdbqt \
    --ligand data/ligands/flumazenil.sdf \
    --output test_out.sdf \
    --profile thorough \
    --center-x 0.0 \
    --center-y 0.0 \
    --center-z 0.0 \
    --size-x 25.0 \
    --size-y 25.0 \
    --size-z 25.0
```

### Validation Run
```bash
python scripts/pipeline/run_validation.py \
    --receptor 6X3U_apo \
    --ligand flumazenil \
    --native data/ligands/flumazenil.sdf \
    --config config.json \
    --out validation_report.md
```

### Production Dispatch
```bash
python scripts/pipeline/local_dispatcher.py \
    --config scripts/pipeline/local_dispatcher.json \
    --gnina-profile thorough \
    --creds credentials/key.json
```

## Architecture Changes

### Old (Broken)
- GNINA parameters hardcoded in Kaggle worker
- exhaustiveness=8, num_modes=9, cnn_scoring=rescore
- No validation workflow
- No score quality assessment

### New (Fixed)
- Centralized configuration in `gnina_runner.py`
- Profile-based parameter selection
- Automated validation pipeline
- Quality flags in score extraction
- Reusable worker architecture

## File Structure

```
scripts/pipeline/
  gnina_runner.py           - Core GNINA execution logic
  run_validation.py         - Automated validation workflow
  extract_scores.py         - Score extraction with quality checks
  local_dispatcher.py       - Job queue builder for Kaggle

scripts/kaggle/
  worker.py                 - Generic multi-GPU task worker
  gnina_worker_deprecated.py - Old hardcoded version (backup)
```

## Migration Guide

### For Existing Projects

1. Update `work_queue.json` format:
```json
{
  "task_type": "gnina_docking",
  "gnina_profile": "thorough",
  "max_workers": 2,
  "jobs": [
    { "receptor": "...", "ligand": "...", ... }
  ]
}
```

2. Update Kaggle kernel to use new worker:
```python
# In kernel-metadata.json, ensure worker.py is listed
# Old: gnina_worker.py
# New: worker.py
```

3. Ensure pipeline module is accessible:
```python
# In worker.py
sys.path.append("/kaggle/working/scripts")
from pipeline.gnina_runner import GNINARunner, GNINAConfig
```

## Recommendations

### Pre-Production Checklist
- [ ] Run validation on at least one test case
- [ ] Verify RMSD < 2.0A
- [ ] Confirm negative binding energies
- [ ] Check runtime > 5 minutes
- [ ] Review score quality flags in CSV

### Production Settings
- Use `thorough` profile (exhaustiveness=32)
- Enable quality assessment in score extraction
- Monitor for positive energy warnings
- Save validation reports for documentation

### Future Improvements
- Custom CNN model training for specific targets
- Multi-ligand batching to reduce receptor reloading
- Automated parameter optimization based on validation results
- Integration with molecular dynamics for refinement
