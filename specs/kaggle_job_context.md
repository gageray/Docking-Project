# Kaggle Docking Job - Current Status

## What Was Fixed

### 1. Authentication System
- **build_kaggle_push.py** now uses OAuth from config.json (like rest of pipeline)
- Removed service_account authentication requirement
- Uses drive_auth.py's setup_drive() method

### 2. Work Queue Format
- Updated to include metadata: `task_type`, `gnina_profile`, `max_workers`
- Jobs wrapped in `"jobs"` array
- Profile set to "validation" (exhaustiveness=32, num_modes=20, cnn_scoring=all)

### 3. Script Upload Structure
- Creates `scripts/` folder in Drive push folder
- Uploads to flat structure (no nested pipeline/ folder):
  - `scripts/setup_env.py`
  - `scripts/drive_auth.py`
  - `scripts/drive_io.py`
  - `scripts/worker.py` (NEW - not gnina_worker.py)
  - `scripts/gnina_runner.py` (from pipeline/)

### 4. Worker Import Fix
- **worker.py** line 18: Changed from `from pipeline.gnina_runner import` to `from gnina_runner import`
- Matches flat scripts/ folder structure

### 5. Bootstrap Updates
- Removed hardcoded DRIVE_PUSH_FOLDER and DRIVE_OUTPUTS_FOLDER constants
- Now uses drive_auth.DEFAULT_FOLDERS["push"] and ["outputs"]
- Imports `worker` not `gnina_worker`

### 6. Drive Auth Updates
- Added "push" folder to DEFAULT_FOLDERS dict
- `"push": "1_rJoAgziylRUQhWqXB2deTU1YW5_hBXM"`

---

## Current Job Configuration

**File:** `flumazenil_6x3u_job.json`
```json
{
  "receptors": ["6X3U_apo.pdbqt"],
  "ligands": ["flumazenil.sdf"],
  "boxes": [
    {
      "center": [0.0, 0.0, 0.0],
      "size": [25.0, 25.0, 25.0]
    }
  ]
}
```

---

## Drive Push Folder Contents (VERIFIED CORRECT)

```
push/ (1_rJoAgziylRUQhWqXB2deTU1YW5_hBXM)
├── receptors/
│   └── 6X3U_apo.pdbqt
├── ligands/
│   └── flumazenil.sdf
├── scripts/
│   ├── setup_env.py
│   ├── drive_auth.py
│   ├── drive_io.py
│   ├── worker.py (NEW)
│   └── gnina_runner.py
└── work_queue.json
```

**work_queue.json contents:**
```json
{
  "task_type": "gnina_docking",
  "gnina_profile": "validation",
  "max_workers": 2,
  "jobs": [
    {
      "receptor": "6X3U_apo.pdbqt",
      "ligand": "flumazenil.sdf",
      "cx": 0.0, "cy": 0.0, "cz": 0.0,
      "sx": 25.0, "sy": 25.0, "sz": 25.0
    }
  ]
}
```

---

## Workflow Commands

### 1. Upload Job to Drive
```bash
python build_kaggle_push.py --config flumazenil_6x3u_job.json
```

### 2. Push Kernel to Kaggle
```bash
cd scripts/kaggle
kaggle kernels push
```

**Note:** Only bootstrap.py is pushed to Kaggle. Everything else downloads from Drive push folder.

---

## Problem Encountered

### Issue
Kaggle ran version 7 with OLD code showing:
- Wrong parameters: `--exhaustiveness 8 --num_modes 9 --cnn_scoring rescore`
- 2 jobs instead of 1 (flumazenil + diazepam)
- Old output format: "GNINA Docking Worker" not "Multi-GPU Task Worker"

### But Drive Has Correct Files
- work_queue.json: 1 job, validation profile ✓
- scripts/worker.py: NEW worker (6793 bytes) ✓
- scripts/gnina_runner.py: Present ✓

### Conclusion
Kaggle kernel version 7 ran BEFORE changes were made. Version was already queued/running.

---

## Next Steps

1. **Check Kaggle documentation** - How does kernel versioning work?
   - Does kaggle kernels push create a new version or update existing?
   - How to force a fresh run vs cached session?
   - How to cancel running kernel?

2. **Verify kernel metadata** - Check scripts/kaggle/kernel-metadata.json
   ```json
   {
     "id": "ineptrobot/gnina-docking-job",
     "code_file": "bootstrap.py",
     "kernel_type": "script",
     "enable_gpu": "true"
   }
   ```

3. **Trigger new run** - Need to either:
   - Push version 8 and wait for it to run
   - Cancel current run and restart
   - Check if auto-run is enabled

---

## Key Files Modified

1. **build_kaggle_push.py**
   - Line 17: Import drive_auth
   - Line 21: Use DEFAULT_FOLDERS["push"]
   - Line 24-27: setup_drive() uses OAuth
   - Line 121-126: work_queue with metadata
   - Line 136-157: Upload scripts/ with gnina_runner.py

2. **scripts/kaggle/worker.py**
   - Line 18: `from gnina_runner import` (not pipeline.gnina_runner)

3. **scripts/kaggle/bootstrap.py**
   - Line 19: Use DEFAULT_FOLDERS["push"]
   - Line 74: `import worker` (not gnina_worker)
   - Line 87: Use DEFAULT_FOLDERS["outputs"]
   - Line 140: Use DEFAULT_FOLDERS["outputs"]

4. **scripts/kaggle/drive_auth.py**
   - Line 16: Added "push" folder to DEFAULT_FOLDERS

5. **scripts/pipeline/local_dispatcher.py**
   - Line 113: Upload to DEFAULT_FOLDERS["push"] (not "root")

---

## Expected Behavior (When Fresh Run Happens)

1. Bootstrap downloads from push folder
2. Installs dependencies via setup_env.py
3. Imports worker.py (NEW)
4. worker.py reads work_queue.json metadata
5. Loads GNINAConfig.from_profile("validation")
6. Runs with correct parameters:
   - exhaustiveness=32
   - num_modes=20
   - cnn_scoring=all
7. Processes ONLY flumazenil (1 job)
8. Output: 6X3U_apo_flumazenil_out.sdf
9. Uploads to Drive outputs/ folder

---

## URL
https://www.kaggle.com/code/ineptrobot/gnina-docking-job

Current version: 7 (OLD CODE)
Need to trigger: version 8+ (NEW CODE)
