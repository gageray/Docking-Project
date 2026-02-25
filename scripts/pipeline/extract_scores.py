"""
Demo script to extract GNINA / Vina energies directly from an .sdf file.
Run this locally on your PC after downloading the outputs from Google Drive!
"""
from rdkit import Chem
from pathlib import Path
import csv
import sys

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
            
            results.append({
                "File": Path(sdf_path).name,
                "Pose": pose_num,
                "Vina_Affinity_kcal_mol": props.get("minimizedAffinity", "N/A"),
                "CNN_Score_Probability": props.get("CNNscore", "N/A"),
                "CNN_Affinity_kcal_mol": props.get("CNNaffinity", "N/A"),
                "CNN_Variance": props.get("CNNvariance", "N/A")
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
