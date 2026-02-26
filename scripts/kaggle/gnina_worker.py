import subprocess
import concurrent.futures
import json
import os
import sys
from queue import Queue
from typing import Dict, Any

# Ensure we can import the local Kaggle scripts after cloning
sys.path.insert(0, '/kaggle/working/repo/scripts/kaggle')

# Global queue for thread-safe GPU ID distribution allocation
job_queue: Queue[int] = Queue()

def process_job(job_data: Dict[str, Any]) -> None:
    """
    Given a resolved JSON dictionary, run GNINA on the isolated GPU
    and immediately upload the result to Drive.

    TODO: Benchmark receptor loading overhead. Current architecture loads
    receptor for every ligand. If loading is significant (grid computation,
    CNN model setup), refactor to multi-ligand format where one GNINA call
    handles multiple ligands against same receptor (reduces reloading).
    Test: Time multiple single-ligand runs vs one multi-ligand run.
    """
    global job_queue
    gpu_id = job_queue.get()
    try:
        # Isolated execution
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        
        # Pull parameters from the pre-resolved dictionary
        receptor = f"/kaggle/working/receptors/{job_data['receptor']}"
        ligand = f"/kaggle/working/ligands/{job_data['ligand']}"
        
        # Build the exact output filename
        out_name = f"{job_data['receptor'].replace('.pdbqt', '')}_{job_data['ligand'].replace('.sdf', '')}_out.sdf"
        out_path = f"/kaggle/working/outputs/{out_name}"
        
        print(f"[*] Thread {gpu_id} starting GNINA for {out_name}")
        
        # Ensure the file doesn't already exist (idempotency check per thread)
        if os.path.exists(out_path):
             print(f"[*] Output {out_name} already exists. Skipping GNINA run.")
        else:
            # Run standard GNINA
            cmd = [
                "gnina",
                "-r", receptor,
                "-l", ligand,
                "-o", out_path,
                "--center_x", str(job_data['cx']),
                "--center_y", str(job_data['cy']),
                "--center_z", str(job_data['cz']),
                "--size_x", str(job_data['sx']),
                "--size_y", str(job_data['sy']),
                "--size_z", str(job_data['sz']),
                "--exhaustiveness", "8",
                "--num_modes", "9",
                "--cnn_scoring", "rescore"
            ]
            
            subprocess.run(cmd, env=env, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            
        print(f"[+] Thread {gpu_id} completed {out_name}. Uploading to Drive...")
        
        # As soon as it finishes, upload it to Drive to safely persist it
        import drive_io, drive_auth # type: ignore
        service, _ = drive_auth.setup_drive("/kaggle/input/my-creds/key.json", verify=False)
        drive_io.upload_file(service, out_path, drive_auth.DEFAULT_FOLDERS["outputs"])
        
    except Exception as e:
        print(f"[-] Job failed on GPU {gpu_id}: {e}")
    finally:
        job_queue.put(gpu_id)

def bootstrap_environment():
    """Clone the repo and install dependencies before running jobs."""
    print("=== Bootstrapping Environment ===")
    
    # 1. We must download the clone repo script from GitHub directly using standard tools since it's a chicken-egg problem.
    os.makedirs("/kaggle/working/repo", exist_ok=True)
    
    # Or instead, assume we are pushing the kaggle directory including `gnina_worker.py`, `clone_repo.py`, `setup_env.py` directly as kernel files
    # Yes, the dispatcher pushes `scripts/kaggle/*`! So they are in `sys.path` where the kernel runs.
    # We can just import them directly!
    
    import clone_repo # type: ignore
    print("[*] Cloning repository...")
    repo_info = clone_repo.setup_repo(
        repo_url="https://github.com/gageray/Docking-Project.git",
        dest_path="/kaggle/working/repo",
        pull_if_exists=True
    )
    
    # 2. Add repo Kaggle scripts to path just in case
    sys.path.insert(0, '/kaggle/working/repo/scripts/kaggle')
    
    import setup_env # type: ignore
    print("[*] Installing environment tools...")
    setup_env.setup_kaggle_environment(verify=True)
    
    # 3. Drive Authenticate and pull Data
    print("[*] Contacting Google Drive for data download...")
    import drive_auth, drive_io # type: ignore
    
    # Kaggle dataset secrets MUST be mounted here for this to work
    key_path = "/kaggle/input/my-creds/key.json"
    if not os.path.exists(key_path):
        print(f"[-] WARNING: Drive credentials not found at {key_path}")
        print("[-] Execution will likely fail.")
        
    service, _ = drive_auth.setup_drive(key_path, verify=False)
    
    # Define outputs directory so GNINA doesn't complain
    os.makedirs("/kaggle/working/outputs", exist_ok=True)
    
    print("[*] Downloading receptors...")
    drive_io.download_folder(service, drive_auth.DEFAULT_FOLDERS["receptors"], "/kaggle/working/receptors")
    
    print("[*] Downloading ligands...")
    drive_io.download_folder(service, drive_auth.DEFAULT_FOLDERS["ligands"], "/kaggle/working/ligands")
    
    print("[*] Downloading work_queue.json...")
    # work_queue.json should be in root
    files = drive_io.list_drive_folder(service, drive_auth.DEFAULT_FOLDERS["root"], verbose=False)
    queue_file_id = None
    for f in files:
        if f['name'] == 'work_queue.json':
            queue_file_id = f['id']
            break
            
    if queue_file_id:
        drive_io.download_from_drive(service, queue_file_id, "/kaggle/working/work_queue.json")
    else:
        print("[-] FATAL: work_queue.json not found in Drive root.")
        sys.exit(1)


def main() -> None:
    # Bootstrap the Kaggle environment (Git clone, installations, Drive downloads)
    bootstrap_environment()
    
    work_file = '/kaggle/working/work_queue.json'
    if not os.path.exists(work_file):
        print("[-] work_queue.json not found. Exiting.")
        sys.exit(1)
        
    with open(work_file, 'r') as f:
        jobs = json.load(f)
        
    global job_queue
    job_queue.put(0) # GPU 0
    job_queue.put(1) # GPU 1
    
    print(f"[*] Starting concurrent Executor for {len(jobs)} jobs...")
    
    # 2 Workers, matching the 2 T4 GPUs. This intentionally blocks until empty.
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        executor.map(process_job, jobs)
        
    print("[*] All jobs processed and uploaded successfully.")
    sys.exit(0)
        
if __name__ == "__main__":
    main()
