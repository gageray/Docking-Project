#!/usr/bin/env python3
"""GNINA docking worker - processes work queue on multiple GPUs."""
import subprocess
import concurrent.futures
import json
import os
from queue import Queue
from typing import Dict, Any

# Global queue for thread-safe GPU allocation
job_queue: Queue[int] = Queue()

def process_job(job_data: Dict[str, Any]) -> None:
    """Run GNINA on isolated GPU."""
    global job_queue
    gpu_id = job_queue.get()
    try:
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

        receptor = f"receptors/{job_data['receptor']}"
        ligand = f"ligands/{job_data['ligand']}"
        out_name = f"{job_data['receptor'].replace('.pdbqt', '')}_{job_data['ligand'].replace('.sdf', '')}_out.sdf"
        out_path = f"/kaggle/working/outputs/{out_name}"

        print(f"[*] GPU {gpu_id} starting: {out_name}")

        if os.path.exists(out_path):
            print(f"[*] Output exists, skipping: {out_name}")
        else:
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
            subprocess.run(cmd, env=env, check=True)

        print(f"[+] GPU {gpu_id} completed: {out_name}")

    except Exception as e:
        print(f"[-] Job failed on GPU {gpu_id}: {e}")
    finally:
        job_queue.put(gpu_id)

def main() -> None:
    print("=== GNINA Docking Worker ===")

    # Create output directory
    os.makedirs("/kaggle/working/outputs", exist_ok=True)

    # Load work queue
    if not os.path.exists("work_queue.json"):
        print("[-] work_queue.json not found!")
        return

    with open("work_queue.json", 'r') as f:
        jobs = json.load(f)

    print(f"[*] Processing {len(jobs)} jobs on 2 GPUs...")

    # Initialize GPU queue
    global job_queue
    job_queue.put(0)
    job_queue.put(1)

    # Process jobs in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        executor.map(process_job, jobs)

    print("[+] All jobs completed!")

if __name__ == "__main__":
    main()
