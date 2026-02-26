#!/usr/bin/env python3
"""Build a Kaggle push directory with all necessary files for docking jobs.

Usage:
    python build_kaggle_push.py --receptor 6X3U_apo.pdbqt --ligand flumazenil.sdf --box 0,0,0,20,20,20
    python build_kaggle_push.py --config job_config.json
"""

import json
import shutil
import argparse
from pathlib import Path


def build_push_directory(receptors, ligands, boxes, output_dir="kaggle_push", creds_path="service-account-key.json"):
    """
    Build a Kaggle push directory with everything needed for docking.

    Args:
        receptors: List of receptor filenames (e.g., ["6X3U_apo.pdbqt"])
        ligands: List of ligand filenames (e.g., ["flumazenil.sdf"])
        boxes: List of box dicts with center/size (one per receptor)
        output_dir: Where to build the push directory
        creds_path: Path to Google Drive service account key
    """
    root = Path(__file__).parent
    push_dir = root / output_dir

    # Clean and create push directory
    if push_dir.exists():
        shutil.rmtree(push_dir)
    push_dir.mkdir(parents=True)

    print(f"Building Kaggle push directory: {push_dir}")

    # Create subdirs
    (push_dir / "receptors").mkdir()
    (push_dir / "ligands").mkdir()

    # Copy receptor files
    receptor_dir = root / "data" / "receptors" / "prepped"
    for rec in receptors:
        src = receptor_dir / rec
        if not src.exists():
            raise FileNotFoundError(f"Receptor not found: {src}")
        shutil.copy(src, push_dir / "receptors" / rec)
        print(f"  ✓ Copied receptor: {rec}")

    # Copy ligand files
    ligand_dir = root / "data" / "ligands"
    for lig in ligands:
        src = ligand_dir / lig
        if not src.exists():
            raise FileNotFoundError(f"Ligand not found: {src}")
        shutil.copy(src, push_dir / "ligands" / lig)
        print(f"  ✓ Copied ligand: {lig}")

    # Build work queue
    work_queue = []
    for rec, box in zip(receptors, boxes):
        for lig in ligands:
            work_queue.append({
                "receptor": rec,
                "ligand": lig,
                "cx": box["center"][0],
                "cy": box["center"][1],
                "cz": box["center"][2],
                "sx": box["size"][0],
                "sy": box["size"][1],
                "sz": box["size"][2]
            })

    with open(push_dir / "work_queue.json", 'w') as f:
        json.dump(work_queue, f, indent=2)
    print(f"  ✓ Created work_queue.json ({len(work_queue)} jobs)")

    # Copy Kaggle scripts that will be cloned from GitHub
    kaggle_scripts = root / "scripts" / "kaggle"
    needed_scripts = [
        "setup_env.py",
        "drive_auth.py",
        "drive_io.py"
    ]

    for script in needed_scripts:
        src = kaggle_scripts / script
        if src.exists():
            shutil.copy(src, push_dir / script)
            print(f"  ✓ Copied: {script}")

    # Create kernel metadata (bootstrap.py is the entry point)
    metadata = {
        "id": "ineptrobot/gnina-docking-job",
        "title": "GNINA Docking Job",
        "code_file": "bootstrap.py",
        "language": "python",
        "kernel_type": "script",
        "is_private": "true",
        "enable_gpu": "true",
        "enable_internet": "true",
        "dataset_sources": [],
        "competition_sources": [],
        "kernel_sources": [],
        "model_sources": []
    }

    with open(push_dir / "kernel-metadata.json", 'w') as f:
        json.dump(metadata, f, indent=4)
    print(f"  ✓ Created kernel-metadata.json")

    # Copy credentials if they exist
    creds = root / creds_path
    if creds.exists():
        shutil.copy(creds, push_dir / "key.json")
        print(f"  ✓ Copied credentials: key.json")
    else:
        print(f"  ⚠ Credentials not found at {creds_path} - outputs won't upload to Drive")

    # Create gnina_worker.py (simplified - just docking logic)
    worker_script = """#!/usr/bin/env python3
\"\"\"GNINA docking worker - processes work queue on multiple GPUs.\"\"\"
import subprocess
import concurrent.futures
import json
import os
from queue import Queue
from typing import Dict, Any

# Global queue for thread-safe GPU allocation
job_queue: Queue[int] = Queue()

def process_job(job_data: Dict[str, Any]) -> None:
    \"\"\"Run GNINA on isolated GPU.\"\"\"
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
"""

    with open(push_dir / "gnina_worker.py", 'w') as f:
        f.write(worker_script)
    print(f"  ✓ Created gnina_worker.py")

    # Create bootstrap.py (startup/shutdown handler)
    bootstrap_script = """#!/usr/bin/env python3
\"\"\"Bootstrap script for Kaggle notebooks - handles setup, execution, and cleanup.\"\"\"
import subprocess
import sys
import os
import json
import shutil
from datetime import datetime
from pathlib import Path

REPO_URL = "https://github.com/gageray/Docking-Project.git"
SPARSE_PATH = "kaggle_push"

def run_cmd(cmd, check=True, capture=False):
    \"\"\"Run shell command with error handling.\"\"\"
    result = subprocess.run(
        cmd, shell=True, check=check,
        capture_output=capture, text=True
    )
    return result.stdout if capture else None

def clone_workspace():
    \"\"\"Sparse clone just the kaggle_push folder.\"\"\"
    print("\\n=== Cloning workspace from GitHub ===")

    # Clone with sparse checkout
    run_cmd("git clone --depth 1 --filter=blob:none --sparse " + REPO_URL + " /tmp/repo")
    run_cmd("cd /tmp/repo && git sparse-checkout set " + SPARSE_PATH)

    # Copy kaggle_push contents to /kaggle/working
    src = Path(f"/tmp/repo/{SPARSE_PATH}")
    for item in src.iterdir():
        dest = Path("/kaggle/working") / item.name
        if item.is_dir():
            shutil.copytree(item, dest, dirs_exist_ok=True)
        else:
            shutil.copy2(item, dest)

    print(f"  ✓ Workspace cloned to /kaggle/working")

    # Cleanup
    shutil.rmtree("/tmp/repo")

def setup_environment():
    \"\"\"Install dependencies and verify environment.\"\"\"
    print("\\n=== Setting up environment ===")
    sys.path.insert(0, "/kaggle/working")

    import setup_env
    setup_env.setup_kaggle_environment(verify=True)

def run_worker():
    \"\"\"Execute the docking worker.\"\"\"
    print("\\n=== Running worker ===")
    os.chdir("/kaggle/working")

    import gnina_worker
    gnina_worker.main()

def generate_run_id():
    \"\"\"Generate run ID in format YY-MM-DD-XX (hex increment).\"\"\"
    now = datetime.now()
    date_prefix = now.strftime("%y-%m-%d")

    # Check if credentials exist for Drive check
    if not os.path.exists("/kaggle/working/key.json"):
        return f"{date_prefix}-00"

    # Check Drive for existing runs today
    try:
        sys.path.insert(0, "/kaggle/working")
        import drive_auth, drive_io

        service, _ = drive_auth.setup_drive("/kaggle/working/key.json", verify=False)
        files = drive_io.list_drive_folder(
            service,
            drive_auth.DEFAULT_FOLDERS["outputs"],
            verbose=False
        )

        # Find max hex suffix for today
        max_suffix = -1
        for f in files:
            if f['name'].startswith(date_prefix):
                try:
                    suffix = f['name'].split('-')[-1]
                    max_suffix = max(max_suffix, int(suffix, 16))
                except:
                    pass

        next_suffix = max_suffix + 1
        return f"{date_prefix}-{next_suffix:02x}"

    except Exception as e:
        print(f"  ⚠ Could not check Drive for run ID, using 00: {e}")
        return f"{date_prefix}-00"

def create_run_metadata(run_id):
    \"\"\"Create run metadata file.\"\"\"
    metadata = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "worker": "gnina_worker.py"
    }

    # Copy work queue if exists
    if os.path.exists("/kaggle/working/work_queue.json"):
        with open("/kaggle/working/work_queue.json") as f:
            metadata["work_queue"] = json.load(f)

    # Write metadata
    with open("/kaggle/working/outputs/run_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"  ✓ Created run metadata: {run_id}")

def upload_outputs(run_id, zip_output=False):
    \"\"\"Upload outputs folder to Drive.\"\"\"
    print("\\n=== Uploading outputs to Drive ===")

    if not os.path.exists("/kaggle/working/key.json"):
        print("  ⚠ No credentials found, skipping upload")
        return

    try:
        sys.path.insert(0, "/kaggle/working")
        import drive_auth
        from googleapiclient.http import MediaFileUpload

        service, _ = drive_auth.setup_drive("/kaggle/working/key.json", verify=False)
        outputs_folder_id = drive_auth.DEFAULT_FOLDERS["outputs"]

        # Create run folder on Drive
        folder_meta = {
            "name": run_id,
            "mimeType": "application/vnd.google-apps.folder",
            "parents": [outputs_folder_id]
        }
        run_folder = service.files().create(body=folder_meta, fields="id").execute()
        run_folder_id = run_folder["id"]
        print(f"  ✓ Created Drive folder: {run_id}")

        # Upload all files from outputs/
        output_dir = Path("/kaggle/working/outputs")
        for file_path in output_dir.iterdir():
            if file_path.is_file():
                file_meta = {
                    "name": file_path.name,
                    "parents": [run_folder_id]
                }
                media = MediaFileUpload(str(file_path))
                service.files().create(
                    body=file_meta,
                    media_body=media,
                    fields="id"
                ).execute()
                print(f"  ✓ Uploaded: {file_path.name}")

        print(f"  ✓ All outputs uploaded to Drive: {run_id}")

    except Exception as e:
        print(f"  ✗ Upload failed: {e}")

def main():
    print("=" * 60)
    print("KAGGLE DOCKING BOOTSTRAP")
    print("=" * 60)

    try:
        # Startup
        clone_workspace()
        setup_environment()

        # Generate run ID
        run_id = generate_run_id()
        print(f"\\n=== Run ID: {run_id} ===")

        # Execute worker
        run_worker()

        # Shutdown
        create_run_metadata(run_id)
        upload_outputs(run_id)

        print("\\n" + "=" * 60)
        print("BOOTSTRAP COMPLETE")
        print("=" * 60)

    except Exception as e:
        print(f"\\n✗ Bootstrap failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
"""

    with open(push_dir / "bootstrap.py", 'w') as f:
        f.write(bootstrap_script)
    print(f"  ✓ Created bootstrap.py")

    print(f"\n✓ Push directory ready: {push_dir}")
    print(f"\nTo push to Kaggle:")
    print(f"  kaggle kernels push -p {output_dir} --accelerator NvidiaTeslaT4")


def main():
    parser = argparse.ArgumentParser(description="Build Kaggle push directory for docking jobs")
    parser.add_argument("--receptor", nargs="+", help="Receptor filename(s) from data/receptors/prepped/")
    parser.add_argument("--ligand", nargs="+", help="Ligand filename(s) from data/ligands/")
    parser.add_argument("--box", help="Box params: cx,cy,cz,sx,sy,sz (e.g., 0,0,0,20,20,20)")
    parser.add_argument("--config", type=Path, help="JSON config file with receptors/ligands/boxes")
    parser.add_argument("--output", default="kaggle_push", help="Output directory (default: kaggle_push)")
    parser.add_argument("--creds", default="service-account-key.json", help="Path to Drive credentials")

    args = parser.parse_args()

    if args.config:
        # Load from config file
        with open(args.config) as f:
            config = json.load(f)
        receptors = config["receptors"]
        ligands = config["ligands"]
        boxes = config["boxes"]
    elif args.receptor and args.ligand and args.box:
        # Load from CLI args
        receptors = args.receptor
        ligands = args.ligand
        box_vals = [float(x) for x in args.box.split(",")]
        boxes = [{
            "center": box_vals[0:3],
            "size": box_vals[3:6]
        }] * len(receptors)
    else:
        parser.error("Provide either --config or (--receptor + --ligand + --box)")

    build_push_directory(receptors, ligands, boxes, args.output, args.creds)


if __name__ == "__main__":
    main()
