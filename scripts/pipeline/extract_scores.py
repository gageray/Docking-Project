"""
Extract GNINA / Vina energies directly from .sdf files with quality assessment.
Run this locally on your PC after downloading the outputs from Google Drive!
"""
from rdkit import Chem
from pathlib import Path
import csv
import sys


def classify_binding(vina_aff, cnn_aff):
    """Classify binding quality from scores"""
    if vina_aff == "N/A" and cnn_aff == "N/A":
        return "NO_SCORES"

    # Convert to float if available
    try:
        vina = float(vina_aff) if vina_aff != "N/A" else None
        cnn = float(cnn_aff) if cnn_aff != "N/A" else None
    except (ValueError, TypeError):
        return "INVALID_SCORES"

    # Check for positive energies (bad)
    if vina is not None and vina > 0:
        return "POOR_POSITIVE_ENERGY"
    if cnn is not None and cnn > 0:
        return "POOR_POSITIVE_ENERGY"

    # Use best available score
    if vina is not None and cnn is not None:
        best_score = min(vina, cnn)
    elif vina is not None:
        best_score = vina
    elif cnn is not None:
        best_score = cnn
    else:
        return "NO_SCORES"

    if best_score < -10:
        return "EXCELLENT"
    elif best_score < -8:
        return "GOOD"
    elif best_score < -6:
        return "MODERATE"
    else:
        return "WEAK"


def check_warnings(vina_aff, cnn_aff, cnn_score):
    """Flag suspicious results"""
    warnings = []

    # Check for missing Vina score
    if vina_aff == "N/A" or vina_aff == 0.0:
        warnings.append("VINA_NOT_SAVED")

    # Check for positive CNN energy
    try:
        if cnn_aff != "N/A" and float(cnn_aff) > 0:
            warnings.append("POSITIVE_CNN_ENERGY")
    except (ValueError, TypeError):
        pass

    # Check for positive Vina energy
    try:
        if vina_aff != "N/A" and float(vina_aff) > 0:
            warnings.append("POSITIVE_VINA_ENERGY")
    except (ValueError, TypeError):
        pass

    # Check for low CNN confidence
    try:
        if cnn_score != "N/A" and float(cnn_score) < 0.5:
            warnings.append("LOW_CNN_CONFIDENCE")
    except (ValueError, TypeError):
        pass

    return "|".join(warnings) if warnings else "OK"


def parse_sdf_file(sdf_path):
    """Parses a single GNINA _out.sdf file and returns a list of dictionaries for each pose."""
    results = []
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path))
        for i, mol in enumerate(supplier):
            if mol is None:
                continue

            pose_num = i + 1
            props = mol.GetPropsAsDict()

            vina_aff = props.get("minimizedAffinity", "N/A")
            cnn_score = props.get("CNNscore", "N/A")
            cnn_aff = props.get("CNNaffinity", "N/A")
            cnn_var = props.get("CNNvariance", "N/A")

            results.append({
                "File": Path(sdf_path).name,
                "Pose": pose_num,
                "Vina_Affinity_kcal_mol": vina_aff,
                "CNN_Score_Probability": cnn_score,
                "CNN_Affinity_kcal_mol": cnn_aff,
                "CNN_Variance": cnn_var,
                "Binding_Quality": classify_binding(vina_aff, cnn_aff),
                "Warning": check_warnings(vina_aff, cnn_aff, cnn_score)
            })
    except Exception as e:
        print(f"[-] Error reading {sdf_path}: {e}")

    return results

def process_directory(input_dir, output_csv):
    """Scans a directory for all .sdf files and compiles a master CSV table."""
    p_dir = Path(input_dir)
    if not p_dir.exists() or not p_dir.is_dir():
        print(f"[-] Input directory {input_dir} not found.")
        sys.exit(1)
        
    all_results = []
    
    sdf_files = list(p_dir.rglob("*.sdf"))
    print(f"[*] Found {len(sdf_files)} SDF files in {input_dir}. Extracting scores...")
    
    for sdf_file in sdf_files:
        all_results.extend(parse_sdf_file(sdf_file))
        
    if not all_results:
        print("[-] No poses found to extract.")
        return
        
    keys = all_results[0].keys()
    with open(output_csv, 'w', newline='') as f:
        dict_writer = csv.DictWriter(f, keys)
        dict_writer.writeheader()
        dict_writer.writerows(all_results)
        
    print(f"[+] Master table saved: {output_csv} with {len(all_results)} total rows.")

if __name__ == "__main__":
    OUT_DIR = "data/outputs"
    CSV_FILE = "data/analysis/docking_results.csv"
    
    Path("data/analysis").mkdir(parents=True, exist_ok=True)
    process_directory(OUT_DIR, CSV_FILE)
