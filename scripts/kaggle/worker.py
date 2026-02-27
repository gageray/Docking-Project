#!/usr/bin/env python3
"""
Generic multi-GPU task worker for Kaggle kernels.
Supports multiple task types (GNINA docking, MD, etc.) with configurable behavior.
"""
import subprocess
import concurrent.futures
import json
import os
import sys
from queue import Queue
from pathlib import Path
from typing import Dict, Any, Callable

# Add scripts directory to path for imports
sys.path.append("/kaggle/working/scripts")

from pipeline.gnina_runner import GNINARunner, GNINAConfig


# Global queue for thread-safe GPU ID distribution
job_queue: Queue[int] = Queue()


class WorkerConfig:
    """Load worker behavior from work_queue.json metadata"""

    def __init__(self, data: Dict[str, Any]):
        self.task_type = data.get("task_type", "gnina_docking")
        self.gnina_profile = data.get("gnina_profile", "standard")
        self.max_workers = data.get("max_workers", 2)

    def __repr__(self):
        return f"WorkerConfig(task_type={self.task_type}, gnina_profile={self.gnina_profile}, max_workers={self.max_workers})"


class TaskWorker:
    """Generic multi-GPU task processor"""

    def __init__(self, config: WorkerConfig, gpu_queue: Queue):
        self.config = config
        self.gpu_queue = gpu_queue
        self.task_handlers = self._register_handlers()

    def _register_handlers(self) -> Dict[str, Callable]:
        """Map task types to handler functions"""
        return {
            "gnina_docking": self.handle_gnina_docking,
        }

    def handle_gnina_docking(self, job_data: Dict[str, Any]) -> None:
        """
        Execute GNINA docking task.

        TODO: Benchmark receptor loading overhead. Current architecture loads
        receptor for every ligand. If loading is significant (grid computation,
        CNN model setup), refactor to multi-ligand format where one GNINA call
        handles multiple ligands against same receptor (reduces reloading).
        Test: Time multiple single-ligand runs vs one multi-ligand run.
        """
        gpu_id = self.gpu_queue.get()
        try:
            # Load GNINA config from worker settings
            gnina_config = GNINAConfig.from_profile(self.config.gnina_profile)
            gnina_config.gpu_id = gpu_id
            runner = GNINARunner(gnina_config)

            # Build paths
            receptor_path = Path(f"/kaggle/working/receptors/{job_data['receptor']}")
            ligand_path = Path(f"/kaggle/working/ligands/{job_data['ligand']}")

            # Build output filename
            out_name = f"{job_data['receptor'].replace('.pdbqt', '')}_{job_data['ligand'].replace('.sdf', '')}_out.sdf"
            output_path = Path(f"/kaggle/working/outputs/{out_name}")

            # Build box parameters
            box_params = {
                "center_x": job_data['cx'],
                "center_y": job_data['cy'],
                "center_z": job_data['cz'],
                "size_x": job_data['sx'],
                "size_y": job_data['sy'],
                "size_z": job_data['sz'],
            }

            # Set GPU environment
            env_vars = {"CUDA_VISIBLE_DEVICES": str(gpu_id)}

            print(f"[*] GPU {gpu_id} starting GNINA for {out_name}")

            # Run docking
            result = runner.run_docking(
                receptor_path=receptor_path,
                ligand_path=ligand_path,
                output_path=output_path,
                box_params=box_params,
                env_vars=env_vars
            )

            # Check result
            if not result.success:
                print(f"[-] GPU {gpu_id} GNINA failed for {out_name}: {result.error}")
                return

            # Log scores for monitoring
            if result.scores:
                best_vina = min(result.scores, key=lambda s: s.minimized_affinity if s.minimized_affinity != 0 else 999)
                best_cnn = min(result.scores, key=lambda s: s.cnn_affinity if s.cnn_affinity != 0 else 999)

                print(f"[+] GPU {gpu_id} completed: {out_name}")
                print(f"    Runtime: {result.runtime_seconds:.1f}s")
                print(f"    Poses: {result.num_poses}")
                print(f"    Best Vina: {best_vina.minimized_affinity:.2f} kcal/mol")
                print(f"    Best CNN: {best_cnn.cnn_affinity:.2f} kcal/mol (score={best_cnn.cnn_score:.3f})")

                # Warn about suspicious results
                validation = runner.validate_output(output_path)
                if not validation.is_valid:
                    print(f"    WARNING: Validation issues detected:")
                    for warning in validation.warnings:
                        print(f"      - {warning}")
            else:
                print(f"[+] GPU {gpu_id} completed: {out_name} (no scores extracted)")

        except Exception as e:
            print(f"[-] Job failed on GPU {gpu_id}: {e}")
        finally:
            self.gpu_queue.put(gpu_id)

    def process_job(self, job_data: Dict[str, Any]) -> None:
        """Route job to appropriate handler"""
        task_type = job_data.get("task_type", self.config.task_type)
        handler = self.task_handlers.get(task_type)

        if handler:
            handler(job_data)
        else:
            print(f"[-] Unknown task type: {task_type}")


def main() -> None:
    """Process all jobs from work_queue.json using available GPUs."""
    print("=" * 60)
    print("Multi-GPU Task Worker")
    print("=" * 60)

    # Ensure outputs directory exists
    os.makedirs("/kaggle/working/outputs", exist_ok=True)

    # Load work queue
    work_file = '/kaggle/working/work_queue.json'
    if not os.path.exists(work_file):
        print("[-] work_queue.json not found. Exiting.")
        sys.exit(1)

    with open(work_file, 'r') as f:
        data = json.load(f)

    # Extract metadata and jobs
    worker_config = WorkerConfig(data)
    jobs = data.get("jobs", [])

    print(f"[*] Configuration: {worker_config}")
    print(f"[*] Jobs to process: {len(jobs)}")

    if len(jobs) == 0:
        print("[*] No jobs to process. Exiting.")
        sys.exit(0)

    # Initialize GPU queue
    global job_queue
    for gpu_id in range(worker_config.max_workers):
        job_queue.put(gpu_id)

    print(f"[*] Starting concurrent executor with {worker_config.max_workers} workers...")

    # Process jobs
    worker = TaskWorker(worker_config, job_queue)
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker_config.max_workers) as executor:
        executor.map(worker.process_job, jobs)

    print("[*] All jobs processed successfully.")
    sys.exit(0)


if __name__ == "__main__":
    main()
