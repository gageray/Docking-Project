This plan outlines a **Headless Batch Workflow** for Kaggle. It eliminates the web UI, utilizes your existing modular utilities, and implements a watchdog to protect your GPU quota.

---

# Kaggle Headless Docking: Implementation Plan

## 1. System Architecture

This workflow separates **Infrastructure** (Kaggle-specific setup) from **Logic** (your docking pipeline in Git).

| **Component**     | **Location**   | **Role**                                                        |
| ----------------- | -------------- | --------------------------------------------------------------- |
| **Antigravity**   | Remote Server  | Development, Git push, and Kaggle API trigger.                  |
| **GitHub**        | Remote Repo    | Source of truth for your actual docking logic and modules.      |
| **Kaggle Script** | Cloud Instance | Headless worker that runs the `entrypoint.py` once and dies.    |
| **Google Drive**  | Cloud Storage  | Persistent storage for input ligands and final docking results. |

---

## 2. The Entrypoint (`entrypoint.py`)

This is the only file pushed to Kaggle. It orchestrates your uploaded modules.

Python

```
import os
import sys
import time
import signal
import threading
from pathlib import Path

# Initial Bootstrapping Imports
import setup_env
import clone_repo
import drive_auth
import drive_io
import setup_dirs
import hardware_check

# --- WATCHDOG CONFIG ---
IDLE_TIMEOUT_MINS = 20  # Kill if no new files in 20 mins
CREDS_PATH = "/kaggle/input/creds-dataset/key.json" # Update to your dataset path

class DockingWatchdog:
    def __init__(self, output_dir, timeout_mins):
        self.output_dir = Path(output_dir)
        self.timeout_sec = timeout_mins * 60
        self.last_count = 0
        self.last_change = time.time()
        self.active = True

    def monitor(self):
        print(f"[*] Watchdog active ({IDLE_TIMEOUT_MINS}m timeout)")
        while self.active:
            # Check for progress by counting result files
            current_count = len(list(self.output_dir.glob("**/*")))
            if current_count > self.last_count:
                self.last_count = current_count
                self.last_change = time.time()

            if (time.time() - self.last_change) > self.timeout_sec:
                print("!!! WATCHDOG: Idle timeout. Terminating to save GPU.")
                os.kill(os.getpid(), signal.SIGTERM)
            time.sleep(60)

def main():
    # 1. Environment & Directory Setup
    setup_env.setup_kaggle_environment(verify=True)
    dirs = setup_dirs.setup_kaggle_dirs()
    hw = hardware_check.verify_hardware_setup()

    # 2. Sync Logic from Git
    # This pulls the latest docking-core logic you just pushed from Antigravity
    repo = clone_repo.setup_repo(
        repo_url="https://github.com/your-username/docking-core.git",
        dest_path="/kaggle/working/pipeline",
        pull_if_exists=True
    )
    sys.path.append(str(repo['path']))

    # 3. Start Watchdog
    dog = DockingWatchdog(dirs['outputs'], IDLE_TIMEOUT_MINS)
    threading.Thread(target=dog.monitor, daemon=True).start()

    # 4. Authenticate & Execute
    try:
        service, _ = drive_auth.setup_drive(CREDS_PATH)

        # --- EXECUTE YOUR PIPELINE ---
        # from pipeline.engine import run_docking
        # run_docking(gpu_count=hw['gpus']['count'])

    except Exception as e:
        print(f"[-] Execution Failure: {e}")
    finally:
        # 5. Persistence & Auto-Shutdown
        print("[*] Uploading results...")
        drive_io.upload_folder(service, dirs['outputs'], drive_auth.DEFAULT_FOLDERS['outputs'])
        dog.active = False
        print("[*] Shutdown triggered.")
        sys.exit(0)

if __name__ == "__main__":
    main()
```

---

## 3. Deployment Steps

### Step A: One-Time Secret Setup

1. Upload your Google Service Account JSON to a **Private Kaggle Dataset** (e.g., `my-creds`).

2. Note the mount path: `/kaggle/input/my-creds/key.json`.

### Step B: The Metadata File (`kernel-metadata.json`)

Create this in your local project root on Antigravity.

JSON

```
{
  "id": "your-username/docking-worker-gpu",
  "title": "Docking Worker GPU",
  "code_file": "entrypoint.py",
  "language": "python",
  "kernel_type": "script",
  "is_private": "true",
  "enable_gpu": "true",
  "enable_internet": "true",
  "dataset_sources": ["your-username/my-creds"],
  "competition_sources": [],
  "kernel_sources": []
}
```

### Step C: The Development Loop

Run these commands from your Antigravity terminal:

1. **Update Logic:** `git commit -am "tuning params" && git push`

2. **Trigger Worker:** `kaggle kernels push -p .`

3. **Monitor Progress:** `kaggle kernels status your-username/docking-worker-gpu`

4. **Retrieve Logs:** `kaggle kernels output your-username/docking-worker-gpu`

---

## 4. Logic Summary

- **Safety:** The watchdog monitors the filesystem. If GNINA (or your engine) hangs or fails to produce output files within 20 minutes, the script kills itself.

- **Efficiency:** Using `kernel_type: script` ensures that as soon as `sys.exit(0)` is hit, the VM is released. No idle web sessions.

- **Cleanliness:** Your `entrypoint.py` and modular `.py` files handle the "plumbing" (Env, Git, Drive). Your actual research code lives in its own repo and is imported dynamically.

**Next Step:** Would you like a single bash script for Antigravity that automates the `git push` and `kaggle push` sequence into one command?
