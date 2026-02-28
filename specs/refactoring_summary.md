# Pipeline Refactoring Summary

## What Was Done

### Files Created
1. **scripts/pipeline/gnina_runner.py** - GNINA execution module with configuration profiles
2. **scripts/pipeline/run_validation.py** - Automated validation workflow
3. **specs/gnina_configuration.md** - GNINA parameter documentation

### Files Renamed
1. **scripts/kaggle/gnina_worker.py → scripts/kaggle/worker.py**
   - Old version backed up as `gnina_worker_deprecated.py`
   - Original `gnina_worker.py` deleted

### Files Modified
1. **scripts/pipeline/local_dispatcher.py** - Added metadata support
2. **scripts/pipeline/extract_scores.py** - Added quality assessment
3. **scripts/pipeline/local_dispatcher.json** - Added new config parameters
4. **config.json** - Added validation and production sections
5. **specs/fixing_docking_validation.md** - Added resolution section

---

## How The System Works

### Architecture Overview

```
LOCAL MACHINE                           GOOGLE DRIVE                    KAGGLE
-------------                           ------------                    ------
data/receptors/prepped/
  └─ 6X3U_apo.pdbqt
data/ligands/
  └─ flumazenil.sdf

        │
        │ local_dispatcher.py scans files
        ├─ Builds work_queue.json
        │  {
        │    "task_type": "gnina_docking",
        │    "gnina_profile": "thorough",
        │    "max_workers": 2,
        │    "jobs": [
        │      {"receptor": "6X3U_apo.pdbqt",
        │       "ligand": "flumazenil.sdf",
        │       "cx": 0.0, "cy": 0.0, "cz": 0.0,
        │       "sx": 25.0, "sy": 25.0, "sz": 25.0}
        │    ]
        │  }
        │
        ├─ Uploads work_queue.json ────────────> Drive: work_queue.json
        │
        └─ Pushes Kaggle kernel ────────────────────────────────────────> Kernel starts
                                                                               │
                                                                               │ bootstrap.py
                                                                               ├─ Downloads work_queue.json
                                                                               ├─ Downloads receptors/
                                                                               ├─ Downloads ligands/
                                                                               │
                                                                               │ worker.py runs
                                                                               ├─ Imports gnina_runner.py
                                                                               ├─ Loads GNINAConfig.from_profile("thorough")
                                                                               ├─ Processes jobs with 2 GPU workers
                                                                               │
                                                                               │ GNINA executes
                                                                               ├─ Receptor: /kaggle/working/receptors/6X3U_apo.pdbqt
                                                                               ├─ Ligand: /kaggle/working/ligands/flumazenil.sdf
                                                                               ├─ Output: /kaggle/working/outputs/6X3U_apo_flumazenil_out.sdf
                                                                               │
                                                                               └─ Uploads outputs/ ──────> Drive: outputs/
```

### Job Dispatch Flow

**Step 1: Local Setup**
- Receptors in `data/receptors/prepped/*.pdbqt`
- Ligands in `data/ligands/*.sdf`
- Box parameters in `box_params.json` (or use defaults)

**Step 2: Dispatch**
```bash
python scripts/pipeline/local_dispatcher.py \
    --config scripts/pipeline/local_dispatcher.json \
    --gnina-profile thorough
```

What happens:
1. Scans `data/receptors/prepped/` for `*_apo.pdbqt` files
2. Scans `data/ligands/` for `*.sdf` files
3. Cross-products all receptor-ligand pairs
4. Polls Google Drive to skip already-completed jobs
5. Builds `work_queue.json` with metadata + job list
6. Uploads `work_queue.json` to Drive
7. Pushes Kaggle kernel via `kaggle kernels push -p scripts/kaggle`

**Step 3: Kaggle Execution**
- Kernel starts on Kaggle with 2x T4 GPUs
- `bootstrap.py` sets up environment
- Downloads `work_queue.json` from Drive
- Downloads receptors and ligands from Drive
- `worker.py` runs:
  - Loads metadata: task_type, gnina_profile, max_workers
  - Imports `gnina_runner.py` from `/kaggle/working/scripts/pipeline/`
  - Loads `GNINAConfig.from_profile("thorough")` → exhaustiveness=32, num_modes=20, cnn_scoring=all
  - Spawns 2 worker threads (one per GPU)
  - Each thread:
    - Grabs GPU ID from queue
    - Runs GNINA via GNINARunner.run_docking()
    - Logs scores and validation warnings
    - Puts GPU ID back in queue
- Outputs saved to `/kaggle/working/outputs/`
- Outputs uploaded to Drive

**Step 4: Download Results**
```bash
# Pull outputs from Drive to local
python scripts/kaggle/pull_outputs.py

# Extract scores with quality assessment
python scripts/pipeline/extract_scores.py
```

---

## How GNINA Runner Works

### Module Structure

**gnina_runner.py** contains:

1. **GNINAConfig** - Immutable configuration dataclass
   - 4 profiles: `quick_screening`, `standard`, `thorough`, `validation`
   - Parameters: exhaustiveness, num_modes, cnn_scoring, cpu, gpu_id

2. **PoseScore** - Individual pose scoring data
   - Fields: pose_number, minimized_affinity, cnn_score, cnn_affinity, cnn_vs, cnn_variance
   - Methods: `is_good_binding()`, `has_positive_energy()`

3. **ValidationReport** - Output quality check results
   - Fields: file_exists, is_readable, num_poses, has_scores, has_negative_energies, warnings
   - Property: `is_valid` (overall pass/fail)

4. **DockingResult** - Complete run results
   - Fields: success, output_file, runtime_seconds, num_poses, scores, error

5. **GNINARunner** - Main execution class
   - `run_docking()` - Execute GNINA subprocess
   - `extract_scores()` - Parse SDF output
   - `validate_output()` - Quality checks

### Configuration Profiles

| Profile | Exhaustiveness | Num Modes | CNN Scoring | Use Case |
|---------|----------------|-----------|-------------|----------|
| quick_screening | 8 | 9 | none | Initial filter (1-2 min) |
| standard | 16 | 20 | all | Production (5-10 min) |
| thorough | 32 | 20 | all | Final candidates (10-20 min) |
| validation | 32 | 20 | all | Benchmarking (10-20 min) |

**Critical Fix Applied:**
- OLD: exhaustiveness=8, num_modes=9, cnn_scoring="rescore" (BROKEN)
- NEW: exhaustiveness=32, num_modes=20, cnn_scoring="all" (FIXED)

### Usage in Worker

```python
# worker.py imports gnina_runner
from pipeline.gnina_runner import GNINARunner, GNINAConfig

# Load configuration from work_queue.json metadata
gnina_config = GNINAConfig.from_profile("thorough")  # or "standard", "validation", etc.

# Create runner
runner = GNINARunner(gnina_config)

# Execute docking
result = runner.run_docking(
    receptor_path=Path("/kaggle/working/receptors/6X3U_apo.pdbqt"),
    ligand_path=Path("/kaggle/working/ligands/flumazenil.sdf"),
    output_path=Path("/kaggle/working/outputs/6X3U_apo_flumazenil_out.sdf"),
    box_params={"center_x": 0.0, "center_y": 0.0, "center_z": 0.0,
                "size_x": 25.0, "size_y": 25.0, "size_z": 25.0},
    env_vars={"CUDA_VISIBLE_DEVICES": "0"}
)

# Check results
if result.success:
    print(f"Generated {result.num_poses} poses")
    best = min(result.scores, key=lambda s: s.cnn_affinity)
    print(f"Best: {best.cnn_affinity:.2f} kcal/mol")
```

---

## Worker Refactoring

### Old Architecture (gnina_worker.py)

```python
# HARDCODED GNINA parameters
cmd = [
    "gnina",
    "-r", receptor,
    "-l", ligand,
    "-o", out_path,
    "--exhaustiveness", "8",        # HARDCODED
    "--num_modes", "9",             # HARDCODED
    "--cnn_scoring", "rescore"      # HARDCODED (loses Vina scores)
]
subprocess.run(cmd, ...)
```

**Problems:**
- GNINA parameters hardcoded in worker
- No configuration flexibility
- No validation or error handling
- Only works for GNINA (not reusable)

### New Architecture (worker.py)

```python
# CONFIGURATION-DRIVEN
class WorkerConfig:
    def __init__(self, data: dict):
        self.task_type = data.get("task_type", "gnina_docking")
        self.gnina_profile = data.get("gnina_profile", "standard")
        self.max_workers = data.get("max_workers", 2)

class TaskWorker:
    def handle_gnina_docking(self, job_data):
        # Load config from metadata
        gnina_config = GNINAConfig.from_profile(self.config.gnina_profile)
        runner = GNINARunner(gnina_config)

        # Run docking
        result = runner.run_docking(...)

        # Validate and log
        if not result.success:
            print(f"[-] GNINA failed: {result.error}")

        validation = runner.validate_output(output_path)
        if not validation.is_valid:
            for warning in validation.warnings:
                print(f"  - {warning}")
```

**Improvements:**
- GNINA logic in separate module (gnina_runner.py)
- Configuration from work_queue.json metadata
- Proper error handling and validation
- Extensible to other task types (MD, FEP, etc.)
- Score logging and quality warnings

---

## Validation Pipeline (run_validation.py)

### Purpose
Automated self-docking validation to verify GNINA parameters before production runs.

### Workflow

1. **Run GNINA** with validation-grade parameters
   - Uses `GNINAConfig.validation()` profile
   - exhaustiveness=32, num_modes=20, cnn_scoring=all

2. **Calculate RMSD** vs native pose
   - Uses `score_native_rmsd.py` module
   - Compares docked poses to crystal structure

3. **Extract Scores** from output SDF
   - Uses `GNINARunner.extract_scores()`
   - Parses minimizedAffinity, CNNaffinity, CNNscore, etc.

4. **Generate Report** with pass/fail criteria
   - Markdown format
   - Tables for RMSD and scores
   - Validation criteria checklist
   - Recommendations

### Success Criteria

- [ ] Best RMSD < 2.0 Angstroms
- [ ] Best RMSD pose ranks in top 5 by score
- [ ] At least one pose has negative binding energy
- [ ] Runtime > 300 seconds (sufficient exhaustiveness)

### Usage

```bash
python scripts/pipeline/run_validation.py \
    --receptor 6X3U_apo \
    --ligand flumazenil \
    --native data/ligands/flumazenil.sdf \
    --config config.json \
    --out validation_report.md
```

**NOTE:** This was designed to run locally but requires the `docking` conda environment (rdkit dependency). Currently NOT TESTED because:
1. Script requires rdkit module
2. Test was attempted in wrong environment
3. User doesn't test locally - jobs are dispatched to Kaggle

---

## Score Extraction Enhancement

### Old Version (extract_scores.py)

```python
# Basic score extraction
results.append({
    "File": sdf_path.name,
    "Pose": pose_num,
    "Vina_Affinity_kcal_mol": props.get("minimizedAffinity", "N/A"),
    "CNN_Score_Probability": props.get("CNNscore", "N/A"),
    "CNN_Affinity_kcal_mol": props.get("CNNaffinity", "N/A"),
    "CNN_Variance": props.get("CNNvariance", "N/A")
})
```

### New Version (extract_scores.py)

```python
# Score extraction WITH quality assessment
vina_aff = props.get("minimizedAffinity", "N/A")
cnn_aff = props.get("CNNaffinity", "N/A")
cnn_score = props.get("CNNscore", "N/A")

results.append({
    "File": sdf_path.name,
    "Pose": pose_num,
    "Vina_Affinity_kcal_mol": vina_aff,
    "CNN_Score_Probability": cnn_score,
    "CNN_Affinity_kcal_mol": cnn_aff,
    "CNN_Variance": props.get("CNNvariance", "N/A"),
    "Binding_Quality": classify_binding(vina_aff, cnn_aff),  # NEW
    "Warning": check_warnings(vina_aff, cnn_aff, cnn_score)   # NEW
})
```

### Quality Classification

**classify_binding():**
- EXCELLENT: < -10 kcal/mol
- GOOD: -8 to -10 kcal/mol
- MODERATE: -6 to -8 kcal/mol
- WEAK: > -6 kcal/mol
- POOR_POSITIVE_ENERGY: > 0 kcal/mol (FAILURE)
- NO_SCORES: Missing data
- INVALID_SCORES: Parse error

**check_warnings():**
- VINA_NOT_SAVED: minimizedAffinity = 0.0 or N/A
- POSITIVE_CNN_ENERGY: CNNaffinity > 0
- POSITIVE_VINA_ENERGY: minimizedAffinity > 0
- LOW_CNN_CONFIDENCE: CNNscore < 0.5

### CSV Output Example

```csv
File,Pose,Vina_Affinity_kcal_mol,CNN_Score_Probability,CNN_Affinity_kcal_mol,CNN_Variance,Binding_Quality,Warning
6X3U_apo_flumazenil_out.sdf,1,-9.2,0.95,-9.5,0.12,EXCELLENT,OK
6X3U_apo_flumazenil_out.sdf,2,0.0,0.85,4.04,0.53,POOR_POSITIVE_ENERGY,VINA_NOT_SAVED|POSITIVE_CNN_ENERGY
```

---

## Configuration Updates

### config.json

**Added sections:**

```json
{
  "docking_validation": {
    "enabled": true,
    "test_cases": [
      {
        "name": "6X3U_flumazenil_self_dock",
        "receptor": "6X3U_apo",
        "ligand": "flumazenil",
        "native_pose": "data/ligands/flumazenil.sdf",
        "expected_rmsd_threshold": 2.0,
        "run_before_production": true
      }
    ],
    "gnina_profile": "validation",
    "output_dir": "data/validation",
    "receptor_dir": "data/receptors/prepped",
    "ligand_dir": "data/ligands"
  },

  "docking_production": {
    "gnina_profile": "thorough",
    "receptor_dir": "data/receptors/prepped",
    "ligand_dir": "data/ligands",
    "output_dir": "data/outputs",
    "box_config": "box_params.json"
  }
}
```

### local_dispatcher.json

**Added parameters:**

```json
{
  "box_config": "scripts/scheduler/box_config.json",
  "receptor_dir": "data/receptors/prepped",
  "ligand_dir": "data/ligands",
  "work_queue": "work_queue.json",
  "creds": "service-account-key.json",
  "kaggle_dir": "scripts/kaggle",
  "gnina_profile": "thorough",        // NEW
  "task_type": "gnina_docking",       // NEW
  "max_workers": 2                    // NEW
}
```

### work_queue.json Format

**New structure with metadata:**

```json
{
  "task_type": "gnina_docking",
  "gnina_profile": "thorough",
  "max_workers": 2,
  "generated_at": "2026-02-27T07:30:00",
  "jobs": [
    {
      "receptor": "6X3U_apo.pdbqt",
      "ligand": "flumazenil.sdf",
      "cx": 0.0,
      "cy": 0.0,
      "cz": 0.0,
      "sx": 25.0,
      "sy": 25.0,
      "sz": 25.0
    }
  ]
}
```

**Old structure (flat array):**

```json
[
  {
    "receptor": "6X3U_apo.pdbqt",
    "ligand": "flumazenil.sdf",
    "cx": 0.0,
    "cy": 0.0,
    "cz": 0.0,
    "sx": 25.0,
    "sy": 25.0,
    "sz": 25.0
  }
]
```

---

## Path Structure

### Local Machine

```
project/
├── data/
│   ├── receptors/
│   │   └── prepped/
│   │       └── 6X3U_apo.pdbqt
│   ├── ligands/
│   │   ├── flumazenil.sdf
│   │   ├── diazepam.sdf
│   │   └── ...
│   ├── outputs/               (downloaded from Drive)
│   │   └── 6X3U_apo_flumazenil_out.sdf
│   └── validation/            (local validation runs)
│       └── test_flumazenil_out.sdf
├── scripts/
│   ├── pipeline/
│   │   ├── gnina_runner.py
│   │   ├── run_validation.py
│   │   ├── local_dispatcher.py
│   │   ├── extract_scores.py
│   │   └── score_native_rmsd.py
│   └── kaggle/
│       ├── worker.py
│       ├── bootstrap.py
│       ├── drive_auth.py
│       ├── drive_io.py
│       └── gnina_worker_deprecated.py
├── config.json
├── box_params.json
└── work_queue.json            (generated by dispatcher)
```

### Kaggle Kernel

```
/kaggle/working/
├── scripts/
│   └── pipeline/
│       └── gnina_runner.py    (uploaded with kernel)
├── receptors/
│   └── 6X3U_apo.pdbqt         (downloaded from Drive)
├── ligands/
│   └── flumazenil.sdf         (downloaded from Drive)
├── outputs/
│   └── 6X3U_apo_flumazenil_out.sdf  (created by GNINA)
└── work_queue.json            (downloaded from Drive)
```

---

## Implementation Issues Encountered

### 1. CLI Test Confusion

**Mistake:** Added `--output` CLI argument to gnina_runner.py and tried to run local test
```bash
python scripts/pipeline/gnina_runner.py \
    --receptor data/receptors/prepped/6X3U_apo.pdbqt \
    --ligand data/ligands/flumazenil.sdf \
    --output data/validation/test_flumazenil_out.sdf
```

**Problem:**
- User doesn't test locally - jobs are dispatched to Kaggle
- Script requires rdkit (not available in current shell environment)
- CLI interface was unnecessary - gnina_runner is imported as module, not run directly

**Resolution:** CLI remains for debugging purposes but is NOT the primary interface

### 2. Path Structure Confusion

**Mistake:** Used local paths (`data/receptors/prepped/`) in context of Kaggle execution

**Reality:**
- **Local:** `data/receptors/prepped/6X3U_apo.pdbqt`
- **Kaggle:** `/kaggle/working/receptors/6X3U_apo.pdbqt`

**Resolution:**
- gnina_runner.py takes Path arguments (works anywhere)
- worker.py hardcodes Kaggle paths
- local_dispatcher.py uses local paths to scan files

### 3. Module Import Dependencies

**Issue:** gnina_runner.py imports rdkit (not available in base environment)

**Solution:**
- Module runs in Kaggle kernel (rdkit available)
- For local use: must activate `docking` conda environment
- Worker imports module after environment setup

---

## How to Use the New System

### Dispatch Jobs to Kaggle

```bash
# 1. Ensure receptors and ligands are ready
ls data/receptors/prepped/*.pdbqt
ls data/ligands/*.sdf

# 2. (Optional) Create box_params.json with binding site coordinates
# If not provided, uses default: center=(0,0,0), size=(25,25,25)

# 3. Run dispatcher
python scripts/pipeline/local_dispatcher.py \
    --config scripts/pipeline/local_dispatcher.json \
    --gnina-profile thorough \
    --creds credentials/key.json

# Dispatcher will:
# - Scan for receptors and ligands
# - Build work_queue.json
# - Upload to Drive
# - Push Kaggle kernel
```

### Monitor Kaggle Execution

1. Go to Kaggle.com → Code → Your kernel
2. Watch logs for:
   - `[*] GPU 0 starting GNINA for 6X3U_apo_flumazenil_out.sdf`
   - `[+] GPU 0 completed: 6X3U_apo_flumazenil_out.sdf`
   - Runtime, poses, scores
   - Validation warnings

### Download and Analyze Results

```bash
# 1. Pull outputs from Drive
python scripts/kaggle/pull_outputs.py

# 2. Extract scores with quality assessment
python scripts/pipeline/extract_scores.py

# 3. Review CSV
# Look for:
# - Binding_Quality = EXCELLENT/GOOD/MODERATE
# - Warning = OK (no issues)
# - Avoid POOR_POSITIVE_ENERGY or VINA_NOT_SAVED warnings
```

---

## Key Architectural Decisions

### 1. Separation of GNINA Logic from Worker

**Rationale:**
- Worker is task-agnostic (can handle MD, FEP, etc. in future)
- GNINA configuration in reusable module
- Single source of truth for parameters

**Implementation:**
- `gnina_runner.py` = GNINA execution logic
- `worker.py` = Generic multi-GPU task processor
- Worker imports gnina_runner as module

### 2. Configuration-Driven Behavior

**Rationale:**
- No hardcoded parameters
- Easy to switch between quick/standard/thorough/validation
- Configuration in JSON, not code

**Implementation:**
- Profiles defined in GNINAConfig class
- Profile name in work_queue.json metadata
- Worker loads profile: `GNINAConfig.from_profile("thorough")`

### 3. Quality Assessment in Score Extraction

**Rationale:**
- Catch docking failures automatically
- Flag suspicious results (positive energies, missing scores)
- Provide guidance on binding quality

**Implementation:**
- `classify_binding()` function
- `check_warnings()` function
- New CSV columns: Binding_Quality, Warning

### 4. Validation Pipeline for Benchmarking

**Rationale:**
- Verify GNINA parameters before expensive production runs
- Standard practice in computational chemistry
- Automated pass/fail criteria

**Implementation:**
- `run_validation.py` script
- RMSD calculation vs native pose
- Markdown report generation
- Success criteria: RMSD < 2.0A, negative energies, proper runtime

---

## What Was NOT Implemented

### 1. Automated Validation Before Production

**Planned:** `local_dispatcher.py` would run validation automatically if `run_before_production: true`

**Not Done:** Validation script exists but not integrated into dispatcher workflow

**Reason:** Requires local execution with rdkit environment

### 2. Box Parameter Integration

**Planned:** Validation workflow would calculate box from native ligand position

**Not Done:** Box params still loaded from `box_params.json` or defaults

**Reason:** Box calculation is separate module (box_definition.py)

### 3. Local Testing

**Planned:** gnina_runner.py CLI for quick local tests

**Not Done:** CLI exists but not tested/verified

**Reason:** User doesn't test locally - all jobs run on Kaggle

---

## Testing Status

### Tested
- None (no execution performed per user request)

### Ready for Testing
1. **Dispatch to Kaggle:**
   ```bash
   python scripts/pipeline/local_dispatcher.py \
       --config scripts/pipeline/local_dispatcher.json \
       --gnina-profile thorough
   ```

2. **Score Extraction:**
   ```bash
   python scripts/pipeline/extract_scores.py
   ```

### Not Ready for Testing
1. **Validation Pipeline:**
   - Requires rdkit environment
   - Requires box_params.json
   - Not tested

---

## Critical Parameters Fixed

| Parameter | Old (Broken) | New (Fixed) | Impact |
|-----------|--------------|-------------|--------|
| exhaustiveness | 8 | 32 | 4x more thorough search |
| num_modes | 9 | 20 | More diverse poses |
| cnn_scoring | "rescore" | "all" | Preserves Vina scores + CNN scores |

**Expected Results:**
- Runtime: 60s → 10-20 minutes (NORMAL)
- Energies: +4 kcal/mol → -8 to -11 kcal/mol (NEGATIVE)
- RMSD: 2.5A → <2.0A (PASS)

---

## Next Steps

1. Test dispatcher with single flumazenil job
2. Monitor Kaggle execution
3. Verify output scores are negative
4. Run RMSD validation
5. If validation passes, run full production with all ligands

---

## File Dependency Map

```
local_dispatcher.py
├─ imports: drive_auth, drive_io
├─ reads: local_dispatcher.json, box_params.json
├─ scans: data/receptors/prepped/, data/ligands/
├─ generates: work_queue.json
└─ uploads to: Google Drive

worker.py (runs on Kaggle)
├─ imports: pipeline.gnina_runner
├─ reads: work_queue.json
├─ uses: /kaggle/working/receptors/, /kaggle/working/ligands/
├─ generates: /kaggle/working/outputs/*.sdf
└─ (outputs uploaded by bootstrap.py)

gnina_runner.py
├─ imports: rdkit.Chem, subprocess, Path
├─ provides: GNINAConfig, GNINARunner, PoseScore, DockingResult
└─ no file dependencies (takes paths as arguments)

extract_scores.py
├─ imports: rdkit.Chem, csv
├─ scans: data/outputs/*.sdf
└─ generates: data/analysis/docking_results.csv

run_validation.py
├─ imports: gnina_runner, score_native_rmsd
├─ reads: config.json, box_params.json
├─ uses: data/receptors/prepped/, data/ligands/
├─ generates: data/validation/*.sdf
└─ outputs: validation_report.md
```

---

## Summary

**Problem:** GNINA parameters were hardcoded in Kaggle worker with insufficient exhaustiveness, wrong scoring mode, resulting in docking failures.

**Solution:**
1. Created reusable GNINA execution module with configuration profiles
2. Refactored worker to be task-agnostic and configuration-driven
3. Added quality assessment to score extraction
4. Created validation pipeline for benchmarking
5. Updated configurations to use thorough parameters

**Architecture:**
- Local dispatcher scans files and uploads work_queue.json to Drive
- Kaggle worker downloads jobs and executes using gnina_runner module
- Results uploaded to Drive and analyzed locally

**Status:** Code complete, not tested. Ready for Kaggle dispatch.
