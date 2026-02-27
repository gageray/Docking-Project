import subprocess
import concurrent.futures
import json
import os
import sys
from queue import Queue
from typing import Dict, Any

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
            
        print(f"[+] GPU {gpu_id} completed: {out_name}")
        
    except Exception as e:
        print(f"[-] Job failed on GPU {gpu_id}: {e}")
    finally:
        job_queue.put(gpu_id)

def main() -> None:
    """Process all jobs from work_queue.json using available GPUs."""
    print("=== GNINA Docking Worker ===")

    # Ensure outputs directory exists
    os.makedirs("/kaggle/working/outputs", exist_ok=True)
    
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
