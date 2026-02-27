#!/usr/bin/env python3
"""
GNINA execution module - centralized configuration and runner logic.
Provides reusable GNINA docking functionality for both local and distributed execution.
"""
import subprocess
import json
from pathlib import Path
from dataclasses import dataclass
from typing import Dict, List, Optional, Any
from rdkit import Chem
import time


@dataclass
class GNINAConfig:
    """Immutable configuration for GNINA docking runs"""
    exhaustiveness: int
    num_modes: int
    cnn_scoring: str  # "none", "rescore", or "all"
    cpu: int = 8
    gpu_id: Optional[int] = None

    @classmethod
    def quick_screening(cls) -> 'GNINAConfig':
        """Fast screening - minimal parameters for initial filtering"""
        return cls(
            exhaustiveness=8,
            num_modes=9,
            cnn_scoring="none",
            cpu=8
        )

    @classmethod
    def standard(cls) -> 'GNINAConfig':
        """Standard production runs - balanced accuracy/speed"""
        return cls(
            exhaustiveness=16,
            num_modes=20,
            cnn_scoring="all",
            cpu=8
        )

    @classmethod
    def thorough(cls) -> 'GNINAConfig':
        """Thorough search - high accuracy for final candidates"""
        return cls(
            exhaustiveness=32,
            num_modes=20,
            cnn_scoring="all",
            cpu=8
        )

    @classmethod
    def validation(cls) -> 'GNINAConfig':
        """Validation-grade parameters for benchmarking"""
        return cls(
            exhaustiveness=32,
            num_modes=20,
            cnn_scoring="all",
            cpu=8
        )

    @classmethod
    def from_profile(cls, profile_name: str) -> 'GNINAConfig':
        """Load configuration from profile name"""
        profiles = {
            "quick_screening": cls.quick_screening,
            "standard": cls.standard,
            "thorough": cls.thorough,
            "validation": cls.validation,
        }

        if profile_name not in profiles:
            raise ValueError(f"Unknown profile: {profile_name}. Valid: {list(profiles.keys())}")

        return profiles[profile_name]()


@dataclass
class PoseScore:
    """Individual pose scoring results"""
    pose_number: int
    minimized_affinity: float  # Vina score (kcal/mol)
    cnn_score: float  # CNN probability (0-1)
    cnn_affinity: float  # CNN affinity (kcal/mol)
    cnn_vs: Optional[float] = None  # CNN virtual screening score
    cnn_variance: Optional[float] = None  # Uncertainty

    def is_good_binding(self) -> bool:
        """Check if this pose indicates favorable binding"""
        return (self.minimized_affinity < -6.0 or self.cnn_affinity < -6.0)

    def has_positive_energy(self) -> bool:
        """Check for suspicious positive energies"""
        return (self.minimized_affinity > 0 or self.cnn_affinity > 0)


@dataclass
class ValidationReport:
    """Validation results for output quality check"""
    file_exists: bool
    is_readable: bool
    num_poses: int
    has_scores: bool
    has_negative_energies: bool
    warnings: List[str]

    @property
    def is_valid(self) -> bool:
        """Overall validation status"""
        return (self.file_exists and
                self.is_readable and
                self.num_poses > 0 and
                self.has_scores and
                self.has_negative_energies)


@dataclass
class DockingResult:
    """Complete docking run results"""
    success: bool
    output_file: Optional[Path]
    runtime_seconds: float
    num_poses: int
    scores: List[PoseScore]
    error: Optional[str] = None


class GNINARunner:
    """Execute GNINA docking with proper error handling"""

    def __init__(self, config: GNINAConfig):
        self.config = config

    def run_docking(
        self,
        receptor_path: Path,
        ligand_path: Path,
        output_path: Path,
        box_params: Dict[str, float],
        env_vars: Optional[Dict[str, str]] = None
    ) -> DockingResult:
        """
        Run GNINA and return structured results.

        Args:
            receptor_path: Path to receptor PDBQT file
            ligand_path: Path to ligand SDF file
            output_path: Path for output SDF file
            box_params: Dictionary with center_x, center_y, center_z, size_x, size_y, size_z
            env_vars: Optional environment variables (e.g., CUDA_VISIBLE_DEVICES)

        Returns:
            DockingResult with success status, scores, and metadata
        """
        start_time = time.time()

        # Validate inputs
        if not receptor_path.exists():
            return DockingResult(
                success=False,
                output_file=None,
                runtime_seconds=0,
                num_poses=0,
                scores=[],
                error=f"Receptor file not found: {receptor_path}"
            )

        if not ligand_path.exists():
            return DockingResult(
                success=False,
                output_file=None,
                runtime_seconds=0,
                num_poses=0,
                scores=[],
                error=f"Ligand file not found: {ligand_path}"
            )

        # Skip if output already exists
        if output_path.exists():
            runtime = time.time() - start_time
            scores = self.extract_scores(output_path)
            return DockingResult(
                success=True,
                output_file=output_path,
                runtime_seconds=runtime,
                num_poses=len(scores),
                scores=scores,
                error=None
            )

        # Build GNINA command
        cmd = [
            "gnina",
            "-r", str(receptor_path),
            "-l", str(ligand_path),
            "-o", str(output_path),
            "--center_x", str(box_params.get("center_x", 0.0)),
            "--center_y", str(box_params.get("center_y", 0.0)),
            "--center_z", str(box_params.get("center_z", 0.0)),
            "--size_x", str(box_params.get("size_x", 25.0)),
            "--size_y", str(box_params.get("size_y", 25.0)),
            "--size_z", str(box_params.get("size_z", 25.0)),
            "--exhaustiveness", str(self.config.exhaustiveness),
            "--num_modes", str(self.config.num_modes),
            "--cnn_scoring", self.config.cnn_scoring,
            "--cpu", str(self.config.cpu)
        ]

        # Set up environment
        env = env_vars.copy() if env_vars else {}

        try:
            # Run GNINA
            result = subprocess.run(
                cmd,
                env={**subprocess.os.environ, **env},
                check=True,
                capture_output=True,
                text=True
            )

            runtime = time.time() - start_time

            # Extract scores from output
            scores = self.extract_scores(output_path)

            return DockingResult(
                success=True,
                output_file=output_path,
                runtime_seconds=runtime,
                num_poses=len(scores),
                scores=scores,
                error=None
            )

        except subprocess.CalledProcessError as e:
            runtime = time.time() - start_time
            return DockingResult(
                success=False,
                output_file=None,
                runtime_seconds=runtime,
                num_poses=0,
                scores=[],
                error=f"GNINA failed: {e.stderr}"
            )
        except Exception as e:
            runtime = time.time() - start_time
            return DockingResult(
                success=False,
                output_file=None,
                runtime_seconds=runtime,
                num_poses=0,
                scores=[],
                error=f"Unexpected error: {str(e)}"
            )

    def extract_scores(self, sdf_path: Path) -> List[PoseScore]:
        """
        Extract all scores from GNINA output SDF.

        Args:
            sdf_path: Path to GNINA output SDF file

        Returns:
            List of PoseScore objects
        """
        scores = []

        try:
            supplier = Chem.SDMolSupplier(str(sdf_path))

            for i, mol in enumerate(supplier):
                if mol is None:
                    continue

                props = mol.GetPropsAsDict()

                # Extract scores with fallback to 0.0
                minimized_aff = float(props.get("minimizedAffinity", 0.0))
                cnn_score = float(props.get("CNNscore", 0.0))
                cnn_aff = float(props.get("CNNaffinity", 0.0))
                cnn_vs = props.get("CNN_VS")
                cnn_var = props.get("CNNvariance")

                scores.append(PoseScore(
                    pose_number=i + 1,
                    minimized_affinity=minimized_aff,
                    cnn_score=cnn_score,
                    cnn_affinity=cnn_aff,
                    cnn_vs=float(cnn_vs) if cnn_vs is not None else None,
                    cnn_variance=float(cnn_var) if cnn_var is not None else None
                ))

        except Exception as e:
            print(f"[-] Error extracting scores from {sdf_path}: {e}")

        return scores

    def validate_output(self, sdf_path: Path) -> ValidationReport:
        """
        Check if docking succeeded based on output.

        Args:
            sdf_path: Path to GNINA output SDF file

        Returns:
            ValidationReport with quality checks
        """
        warnings = []

        # Check file exists
        if not sdf_path.exists():
            return ValidationReport(
                file_exists=False,
                is_readable=False,
                num_poses=0,
                has_scores=False,
                has_negative_energies=False,
                warnings=["Output file does not exist"]
            )

        # Try to read file
        try:
            scores = self.extract_scores(sdf_path)
        except Exception as e:
            return ValidationReport(
                file_exists=True,
                is_readable=False,
                num_poses=0,
                has_scores=False,
                has_negative_energies=False,
                warnings=[f"Cannot read SDF file: {e}"]
            )

        # Check for poses
        if len(scores) == 0:
            warnings.append("No poses found in output")

        # Check for scores
        has_scores = any(s.minimized_affinity != 0.0 or s.cnn_affinity != 0.0 for s in scores)
        if not has_scores:
            warnings.append("All scores are zero")

        # Check for negative energies
        has_negative = any(s.is_good_binding() for s in scores)
        if not has_negative:
            warnings.append("No poses with negative binding energies")

        # Check for positive energies
        has_positive = any(s.has_positive_energy() for s in scores)
        if has_positive:
            warnings.append("Some poses have positive energies (docking failure)")

        return ValidationReport(
            file_exists=True,
            is_readable=True,
            num_poses=len(scores),
            has_scores=has_scores,
            has_negative_energies=has_negative,
            warnings=warnings
        )


def classify_binding_quality(score: PoseScore) -> str:
    """
    Classify binding quality from scores.

    Args:
        score: PoseScore object

    Returns:
        String classification: EXCELLENT, GOOD, MODERATE, WEAK, POOR_POSITIVE_ENERGY, NO_SCORES
    """
    vina = score.minimized_affinity
    cnn_aff = score.cnn_affinity

    # Check for missing scores
    if vina == 0.0 and cnn_aff == 0.0:
        return "NO_SCORES"

    # Check for positive energies (bad)
    if vina > 0 or cnn_aff > 0:
        return "POOR_POSITIVE_ENERGY"

    # Use best available score
    best_score = min(vina, cnn_aff) if vina != 0.0 and cnn_aff != 0.0 else (vina if vina != 0.0 else cnn_aff)

    if best_score < -10:
        return "EXCELLENT"
    elif best_score < -8:
        return "GOOD"
    elif best_score < -6:
        return "MODERATE"
    else:
        return "WEAK"


if __name__ == "__main__":
    # Simple CLI test interface
    import argparse

    parser = argparse.ArgumentParser(description="Run GNINA docking")
    parser.add_argument("--receptor", required=True, help="Receptor PDBQT file")
    parser.add_argument("--ligand", required=True, help="Ligand SDF file")
    parser.add_argument("--output", required=True, help="Output SDF file")
    parser.add_argument("--profile", default="standard", choices=["quick_screening", "standard", "thorough", "validation"], help="GNINA configuration profile")
    parser.add_argument("--center-x", type=float, required=True, help="Box center X")
    parser.add_argument("--center-y", type=float, required=True, help="Box center Y")
    parser.add_argument("--center-z", type=float, required=True, help="Box center Z")
    parser.add_argument("--size-x", type=float, default=25.0, help="Box size X")
    parser.add_argument("--size-y", type=float, default=25.0, help="Box size Y")
    parser.add_argument("--size-z", type=float, default=25.0, help="Box size Z")
    parser.add_argument("--gpu", type=int, help="GPU ID to use")

    args = parser.parse_args()

    # Load configuration
    config = GNINAConfig.from_profile(args.profile)
    if args.gpu is not None:
        config.gpu_id = args.gpu

    # Run docking
    runner = GNINARunner(config)

    box_params = {
        "center_x": args.center_x,
        "center_y": args.center_y,
        "center_z": args.center_z,
        "size_x": args.size_x,
        "size_y": args.size_y,
        "size_z": args.size_z,
    }

    env_vars = {}
    if args.gpu is not None:
        env_vars["CUDA_VISIBLE_DEVICES"] = str(args.gpu)

    print(f"[*] Running GNINA with profile: {args.profile}")
    print(f"[*] Configuration: exhaustiveness={config.exhaustiveness}, num_modes={config.num_modes}, cnn_scoring={config.cnn_scoring}")

    result = runner.run_docking(
        receptor_path=Path(args.receptor),
        ligand_path=Path(args.ligand),
        output_path=Path(args.output),
        box_params=box_params,
        env_vars=env_vars
    )

    if result.success:
        print(f"[+] Docking completed successfully in {result.runtime_seconds:.1f}s")
        print(f"[+] Generated {result.num_poses} poses")

        if result.scores:
            best = min(result.scores, key=lambda s: s.cnn_affinity if s.cnn_affinity != 0 else s.minimized_affinity)
            print(f"[+] Best pose: #{best.pose_number}")
            print(f"    Vina: {best.minimized_affinity:.2f} kcal/mol")
            print(f"    CNN:  {best.cnn_affinity:.2f} kcal/mol (score={best.cnn_score:.3f})")
            print(f"    Quality: {classify_binding_quality(best)}")
    else:
        print(f"[-] Docking failed: {result.error}")
