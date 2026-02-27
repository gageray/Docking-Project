#!/usr/bin/env python3
"""
Automated self-docking validation workflow.
Runs GNINA with validation-grade parameters, calculates RMSD vs native pose,
and generates validation report with pass/fail criteria.
"""
import argparse
import sys
import json
from pathlib import Path
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional
import math

from gnina_runner import GNINARunner, GNINAConfig, PoseScore
from score_native_rmsd import calculate_ligand_rmsd


@dataclass
class ValidationResult:
    """Complete validation run results"""
    passed: bool
    rmsd_results: List[Tuple[int, float, int]]  # (pose_num, rmsd, n_atoms)
    scores: List[PoseScore]
    report: str
    runtime_seconds: float
    error: Optional[str] = None


class ValidationPipeline:
    """Automated self-docking validation"""

    def __init__(self, config_path: Path):
        self.config = self.load_config(config_path)

    def load_config(self, config_path: Path) -> Dict:
        """Load configuration from JSON file"""
        if not config_path.exists():
            print(f"[-] Config file not found: {config_path}")
            print("[*] Using default configuration")
            return {
                "receptor_dir": Path("data/receptors/prepped"),
                "ligand_dir": Path("data/ligands"),
                "output_dir": Path("data/validation"),
                "box_params": {}
            }

        with open(config_path, 'r') as f:
            data = json.load(f)

        # Extract validation config if exists
        validation_cfg = data.get("docking_validation", {})

        return {
            "receptor_dir": Path(validation_cfg.get("receptor_dir", "data/receptors/prepped")),
            "ligand_dir": Path(validation_cfg.get("ligand_dir", "data/ligands")),
            "output_dir": Path(validation_cfg.get("output_dir", "data/validation")),
            "box_params": validation_cfg.get("box_params", {})
        }

    def load_box_params(self, receptor: str) -> Dict[str, float]:
        """Load box parameters for receptor"""
        # Try to load from box_params.json
        box_file = Path("box_params.json")
        if box_file.exists():
            with open(box_file, 'r') as f:
                data = json.load(f)
                rec_key = f"{receptor}_prepped"
                if rec_key in data:
                    box = data[rec_key]
                    return {
                        "center_x": box["center_x"],
                        "center_y": box["center_y"],
                        "center_z": box["center_z"],
                        "size_x": box["size_x"],
                        "size_y": box["size_y"],
                        "size_z": box["size_z"],
                    }

        # Default box
        return {
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 25.0,
            "size_y": 25.0,
            "size_z": 25.0,
        }

    def run_validation(self, receptor: str, ligand: str, native_pose: Path) -> ValidationResult:
        """
        Complete validation workflow:
        1. Run GNINA with validation-grade parameters
        2. Calculate RMSD for all poses
        3. Extract and validate scores
        4. Generate report
        5. Return pass/fail
        """

        # Ensure output directory exists
        self.config["output_dir"].mkdir(parents=True, exist_ok=True)

        # Step 1: Run docking with thorough settings
        print("[1/4] Running GNINA with validation parameters...")
        gnina_config = GNINAConfig.validation()
        runner = GNINARunner(gnina_config)

        receptor_path = self.config["receptor_dir"] / f"{receptor}.pdbqt"
        ligand_path = self.config["ligand_dir"] / f"{ligand}.sdf"
        output_path = self.config["output_dir"] / f"{receptor}_{ligand}_validation.sdf"

        box_params = self.load_box_params(receptor)

        result = runner.run_docking(
            receptor_path=receptor_path,
            ligand_path=ligand_path,
            output_path=output_path,
            box_params=box_params
        )

        if not result.success:
            return ValidationResult(
                passed=False,
                rmsd_results=[],
                scores=[],
                report="",
                runtime_seconds=result.runtime_seconds,
                error=result.error
            )

        # Step 2: Calculate RMSD
        print("[2/4] Calculating RMSD vs native pose...")
        rmsd_results = calculate_ligand_rmsd(str(native_pose), str(result.output_file))

        # Step 3: Extract scores
        print("[3/4] Extracting docking scores...")
        scores = result.scores

        # Step 4: Generate report
        print("[4/4] Generating validation report...")
        passed = self.evaluate_criteria(rmsd_results, scores, result.runtime_seconds)
        report = self.generate_report(rmsd_results, scores, result, receptor, ligand, passed)

        return ValidationResult(
            passed=passed,
            rmsd_results=rmsd_results,
            scores=scores,
            report=report,
            runtime_seconds=result.runtime_seconds
        )

    def evaluate_criteria(self, rmsd_results: List[Tuple[int, float, int]], scores: List[PoseScore], runtime_seconds: float) -> bool:
        """
        Validation success criteria:
        1. At least 1 pose with RMSD < 2.0A
        2. That pose should rank in top 5 by score
        3. At least 1 pose has NEGATIVE binding energy
        4. Runtime > 300s (indicates sufficient exhaustiveness)
        """

        # Find best RMSD
        valid_rmsd = [(p, r, n) for p, r, n in rmsd_results if not math.isnan(r)]
        if not valid_rmsd:
            print("FAIL: No valid RMSD values")
            return False

        best_rmsd_pose, best_rmsd, _ = min(valid_rmsd, key=lambda x: x[1])

        if best_rmsd >= 2.0:
            print(f"FAIL: Best RMSD {best_rmsd:.2f}A >= 2.0A threshold")
            return False

        # Check if best RMSD pose ranks well
        pose_rank = self.get_pose_rank(best_rmsd_pose, scores)
        if pose_rank > 5:
            print(f"WARNING: Best RMSD pose ranked #{pose_rank} (not in top 5)")
            # Don't fail, but warn

        # Check for negative energies
        has_negative = any(s.is_good_binding() for s in scores)
        if not has_negative:
            print("FAIL: No poses have negative binding energy")
            return False

        # Check runtime
        if runtime_seconds < 300:
            print(f"WARNING: Runtime {runtime_seconds:.1f}s < 300s (may not be thorough enough)")

        print(f"PASS: Best RMSD {best_rmsd:.2f}A, ranked #{pose_rank}")
        return True

    def get_pose_rank(self, pose_num: int, scores: List[PoseScore]) -> int:
        """Get rank of pose by CNN affinity score"""
        # Sort by CNN affinity (or Vina if CNN not available)
        sorted_scores = sorted(scores, key=lambda s: s.cnn_affinity if s.cnn_affinity != 0 else s.minimized_affinity)

        for rank, score in enumerate(sorted_scores, 1):
            if score.pose_number == pose_num:
                return rank

        return 999  # Not found

    def generate_report(self, rmsd_results: List[Tuple[int, float, int]], scores: List[PoseScore], result, receptor: str, ligand: str, passed: bool) -> str:
        """Generate markdown validation report"""

        # Find best RMSD
        valid_rmsd = [(p, r, n) for p, r, n in rmsd_results if not math.isnan(r)]
        if valid_rmsd:
            best_rmsd_pose, best_rmsd, _ = min(valid_rmsd, key=lambda x: x[1])
            pose_rank = self.get_pose_rank(best_rmsd_pose, scores)
        else:
            best_rmsd_pose = 0
            best_rmsd = float('nan')
            pose_rank = 999

        has_negative = any(s.is_good_binding() for s in scores)

        # Format RMSD table
        rmsd_table = ""
        for pose_num, rmsd, n_atoms in rmsd_results[:10]:  # Top 10
            if math.isnan(rmsd):
                quality = "ERROR"
                rmsd_str = "N/A"
            elif rmsd < 1.0:
                quality = "EXCELLENT"
                rmsd_str = f"{rmsd:.3f}"
            elif rmsd < 2.0:
                quality = "GOOD"
                rmsd_str = f"{rmsd:.3f}"
            else:
                quality = "POOR"
                rmsd_str = f"{rmsd:.3f}"

            rmsd_table += f"| {pose_num} | {rmsd_str} | {n_atoms} | {quality} |\n"

        # Format score table
        score_table = ""
        for score in scores[:10]:  # Top 10
            score_table += f"| {score.pose_number} | {score.minimized_affinity:.2f} | {score.cnn_affinity:.2f} | {score.cnn_score:.3f} |\n"

        # Generate recommendations
        recommendations = []
        if not passed:
            if best_rmsd >= 2.0:
                recommendations.append("- RMSD too high. Check box definition covers binding site.")
                recommendations.append("- Verify receptor preparation (protonation, hydrogens).")
                recommendations.append("- Try higher exhaustiveness (64) for more thorough search.")
            if not has_negative:
                recommendations.append("- No negative energies suggests docking failure.")
                recommendations.append("- Check ligand preparation (3D conformer quality).")
                recommendations.append("- Verify PDBQT conversion is correct.")
        else:
            recommendations.append("- Validation passed. Parameters are suitable for production runs.")
            if pose_rank > 5:
                recommendations.append("- Consider training custom CNN scoring for better ranking.")

        recommendations_str = "\n".join(recommendations) if recommendations else "No recommendations."

        report = f"""# Docking Validation Report

## Summary
- **Status**: {'PASS' if passed else 'FAIL'}
- **Receptor**: {receptor}
- **Ligand**: {ligand}
- **Runtime**: {result.runtime_seconds:.1f}s
- **Poses Generated**: {len(scores)}

## RMSD Results (Top 10)
| Pose | RMSD (A) | Atoms | Quality |
|------|----------|-------|---------|
{rmsd_table}

**Best Pose**: #{best_rmsd_pose} with {best_rmsd:.3f}A

## Docking Scores (Top 10)
| Pose | Vina (kcal/mol) | CNN Affinity (kcal/mol) | CNN Score |
|------|-----------------|------------------------|-----------|
{score_table}

## Validation Criteria
- [{'X' if not math.isnan(best_rmsd) and best_rmsd < 2.0 else ' '}] Best RMSD < 2.0A: {best_rmsd:.3f}A
- [{'X' if pose_rank <= 5 else ' '}] Best RMSD pose in top 5: #{pose_rank}
- [{'X' if has_negative else ' '}] Negative binding energies detected
- [{'X' if result.runtime_seconds > 300 else ' '}] Sufficient exhaustiveness (>5min runtime): {result.runtime_seconds:.1f}s

## Recommendations
{recommendations_str}
"""

        return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run docking validation workflow")
    parser.add_argument("--receptor", required=True, help="Receptor name (e.g., 6X3U_apo)")
    parser.add_argument("--ligand", required=True, help="Ligand name (e.g., flumazenil)")
    parser.add_argument("--native", required=True, help="Native pose SDF/PDB")
    parser.add_argument("--config", default="config.json", help="Pipeline config")
    parser.add_argument("--out", help="Output report path")
    args = parser.parse_args()

    pipeline = ValidationPipeline(Path(args.config))
    result = pipeline.run_validation(args.receptor, args.ligand, Path(args.native))

    if result.error:
        print(f"\n[-] Validation failed with error: {result.error}")
        sys.exit(1)

    # Write report
    report_path = args.out or f"validation_{args.receptor}_{args.ligand}.md"
    with open(report_path, 'w') as f:
        f.write(result.report)

    print(f"\n{'='*60}")
    print(f"Validation {'PASSED' if result.passed else 'FAILED'}")
    print(f"Report: {report_path}")
    print(f"{'='*60}")

    sys.exit(0 if result.passed else 1)
