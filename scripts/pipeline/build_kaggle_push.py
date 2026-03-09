#!/usr/bin/env python3
"""Build work queue and push Kaggle kernel.

Reads job config, expands zips, builds work queue from ONLY the specified files,
uploads to Drive, and pushes kernel with specified GPU accelerator.

Usage:
    python build_kaggle_push.py --config flumazenil_6x3u_job.json
    python build_kaggle_push.py --config my_job.json --accelerator NvidiaTeslaP100
"""

import argparse
import subprocess
import sys
import json
import zipfile
import itertools
from pathlib import Path
from datetime import datetime

# Import local_dispatcher only for Drive upload
sys.path.insert(0, str(Path(__file__).parent))
import local_dispatcher


def expand_ligands(ligand_list, ligand_dir):
    """Expand zip files in ligand list to individual SDF files.

    Args:
        ligand_list: List of ligand filenames (may include .zip files)
        ligand_dir: Path to ligands directory

    Returns:
        List of individual SDF filenames
    """
    expanded = []

    for lig in ligand_list:
        lig_path = ligand_dir / lig

        if lig.endswith('.zip') and lig_path.exists():
            # Peek inside zip and get SDF files
            try:
                with zipfile.ZipFile(lig_path, 'r') as zf:
                    sdf_files = [name for name in zf.namelist() if name.endswith('.sdf')]
                    expanded.extend(sdf_files)
                    print(f"      ✓ Found {len(sdf_files)} SDFs in {lig}")
            except zipfile.BadZipFile:
                print(f"      ✗ WARNING: Corrupt zip file {lig}")
        elif lig.endswith('.sdf'):
            # Regular SDF file
            expanded.append(lig)
        else:
            print(f"      ✗ WARNING: Unsupported ligand type: {lig}")

    return expanded


def build_work_queue(receptors, ligands, boxes, gnina_profile, max_workers):
    """Build work queue from specific receptor/ligand lists.

    Args:
        receptors: List of receptor filenames
        ligands: List of ligand filenames (already expanded from zips)
        boxes: List of box dicts (one per receptor)
        gnina_profile: GNINA profile name
        max_workers: Number of GPU workers

    Returns:
        Work queue dict with metadata and jobs
    """
    jobs = []

    for rec, box in zip(receptors, boxes):
        for lig in ligands:
            jobs.append({
                "receptor": rec,
                "ligand": lig,
                "cx": box["center"][0],
                "cy": box["center"][1],
                "cz": box["center"][2],
                "sx": box["size"][0],
                "sy": box["size"][1],
                "sz": box["size"][2]
            })

    return {
        "task_type": "gnina_docking",
        "gnina_profile": gnina_profile,
        "max_workers": max_workers,
        "generated_at": datetime.now().isoformat(),
        "jobs": jobs
    }


def main():
    parser = argparse.ArgumentParser(description="Build work queue and push Kaggle kernel")
    parser.add_argument("--config", type=Path, required=True, help="Job configuration JSON")
    parser.add_argument("--accelerator", type=str, default="NvidiaTeslaT4",
                        help="Kaggle GPU accelerator (NvidiaTeslaT4, NvidiaTeslaP100, etc.)")
    parser.add_argument("--gnina-profile", type=str, help="Override GNINA profile")
    parser.add_argument("--max-workers", type=int, help="Override max GPU workers")

    args = parser.parse_args()

    if not args.config.exists():
        print(f"[-] Config file not found: {args.config}")
        sys.exit(1)

    # Load job config
    with open(args.config) as f:
        config = json.load(f)

    receptors = config["receptors"]
    ligands_raw = config["ligands"]
    boxes = config["boxes"]

    # Get paths
    root = Path(__file__).parent.parent.parent
    ligand_dir = root / "data" / "ligands"
    work_queue_path = Path("work_queue.json")
    creds = config.get("creds", "drive_oauth.json")
    kaggle_dir = Path(config.get("kaggle_dir", "scripts/kaggle"))
    gnina_profile = args.gnina_profile or config.get("gnina_profile", "thorough")
    max_workers = args.max_workers or config.get("max_workers", 2)

    print("=" * 60)
    print("KAGGLE JOB BUILDER & PUSHER")
    print("=" * 60)
    print(f"Config: {args.config}")
    print(f"Accelerator: {args.accelerator}")
    print(f"GNINA Profile: {gnina_profile}")
    print(f"Max Workers: {max_workers}")
    print()

    # Step 1: Expand ligands (peek into zips)
    print("[1/3] Expanding ligands from config...")
    print(f"      Receptors: {receptors}")
    print(f"      Ligands (raw): {ligands_raw}")
    print()

    ligands = expand_ligands(ligands_raw, ligand_dir)

    if not ligands:
        print("[-] No ligands found after expansion!")
        sys.exit(1)

    print(f"      ✓ Expanded to {len(ligands)} ligand(s)")
    print()

    # Step 2: Build work queue
    print("[2/3] Building work queue...")
    job_data = build_work_queue(receptors, ligands, boxes, gnina_profile, max_workers)

    print(f"      ✓ Created {len(job_data['jobs'])} job(s)")
    print(f"        Task type: {job_data['task_type']}")
    print(f"        Profile: {job_data['gnina_profile']}")
    print(f"        Workers: {job_data['max_workers']}")

    # Write work queue locally
    with open(work_queue_path, 'w') as f:
        json.dump(job_data, f, indent=2)
    print(f"      ✓ Wrote {work_queue_path}")
    print()

    # Step 3: Upload work queue to Drive
    print("[3/3] Uploading work queue to Drive...")
    local_dispatcher.upload_queue_to_drive(work_queue_path, creds)
    print("      ✓ Uploaded to Drive push folder")
    print()

    # Step 4: Push kernel to Kaggle with accelerator
    print(f"[4/4] Pushing kernel to Kaggle with {args.accelerator}...")
    try:
        subprocess.run(
            ["kaggle", "kernels", "push", "-p", str(kaggle_dir),
             "--accelerator", args.accelerator],
            check=True
        )
        print("      ✓ Kernel pushed successfully!")
        print()
        print("=" * 60)
        print("✓ JOB SUBMITTED TO KAGGLE")
        print("=" * 60)
        print(f"  Jobs: {len(job_data['jobs'])}")
        print(f"  GPU: {args.accelerator} x {max_workers}")
        print(f"  Profile: {gnina_profile}")
        print()
        print("  Check progress at:")
        print("  https://www.kaggle.com/code/ineptrobot/gnina-docking-job")
        print()
    except subprocess.CalledProcessError as e:
        print(f"[-] Failed to push kernel: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
