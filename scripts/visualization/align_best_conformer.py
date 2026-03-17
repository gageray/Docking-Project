#!/usr/bin/env python3
"""
Extract and align best QM conformer with crystal ligand for PyMOL visualization.
Generates a PyMOL macro (.pml) to load them as easily toggleable objects.
"""

import sys
import json
import zipfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign


def main():
    # Paths
    base_dir = Path(__file__).parent.parent.parent
    results_json = base_dir / "data/analysis/flumazenil_qm_ensemble_rmsd.json"
    pdbqt_zip = base_dir / "data/ligands/flumazenil_qm_ensemble_pdbqt.zip"
    receptor_pdb = base_dir / "data/receptors/6X3U_aligned.pdb"
    output_dir = base_dir / "data/visualization"

    output_dir.mkdir(exist_ok=True)
    qm_output_pdb = output_dir / "best_qm_aligned.pdb"
    pml_script = output_dir / "view_alignment.pml"

    # Load results to find best conformer
    with open(results_json) as f:
        results = json.load(f)

    best_idx = results['results'][0]['conformer_idx']
    best_rmsd = results['results'][0]['rmsd_aligned']
    smiles = results['smiles']

    print(f"Best conformer: #{best_idx} with RMSD={best_rmsd:.3f} Å")

    # ========== 1. Extract and Map Crystal Ligand ==========
    crystal_lines = []
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')) and line[17:20].strip() == 'FYP' and line[21] == 'Z':
                crystal_lines.append(line.rstrip('\n'))

    pdb_block = "\n".join(crystal_lines)
    ref_raw = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
    template = Chem.MolFromSmiles(smiles)
    ref_mol = AllChem.AssignBondOrdersFromTemplate(template, ref_raw)

    # Save isolated crystal ligand for easy PyMOL loading
    crystal_output_pdb = output_dir / "crystal_ligand.pdb"
    Chem.MolToPDBFile(ref_mol, str(crystal_output_pdb))

    # ========== 2. Extract and Map QM Conformer ==========
    with zipfile.ZipFile(pdbqt_zip, 'r') as zf:
        pdbqt_files = sorted([f for f in zf.namelist() if f.endswith('.pdbqt')])
        pdbqt_content = zf.read(pdbqt_files[best_idx]).decode('utf-8')

    clean_pdb = [line[:66] for line in pdbqt_content.splitlines() if line.startswith(('ATOM', 'HETATM'))]
    raw_qm_mol = Chem.MolFromPDBBlock("\n".join(clean_pdb), removeHs=True, sanitize=False)
    qm_mol = AllChem.AssignBondOrdersFromTemplate(template, raw_qm_mol)

    # ========== 3. Align and Save ==========
    # Natively rotates qm_mol coordinates to match ref_mol
    rdMolAlign.AlignMol(qm_mol, ref_mol)
    Chem.MolToPDBFile(qm_mol, str(qm_output_pdb))
    print(f"Saved aligned QM conformer to {qm_output_pdb.name}")

    # ========== 4. Generate PyMOL Macro with Embedded Data ==========
    # Read PDB files
    with open(crystal_output_pdb) as f:
        crystal_pdb_str = f.read()

    with open(qm_output_pdb) as f:
        qm_pdb_str = f.read()

    # Write PML with embedded structures
    with open(pml_script, 'w') as f:
        # Embed crystal ligand
        f.write("cmd.read_pdbstr('''\n")
        f.write(crystal_pdb_str)
        f.write("''', 'crystal_pose')\n\n")

        # Embed QM conformer
        f.write("cmd.read_pdbstr('''\n")
        f.write(qm_pdb_str)
        f.write("''', 'qm_conformer')\n\n")

        # Visualization settings
        f.write("hide everything\n")
        f.write("show sticks, crystal_pose\n")
        f.write("show sticks, qm_conformer\n")
        f.write("color cyan, crystal_pose\n")
        f.write("color magenta, qm_conformer\n")
        f.write("zoom\n")

    print(f"Generated self-contained PyMOL macro: {pml_script.name}")
    print("\nTo visualize, simply run:")
    print(f"pymol {pml_script}")


if __name__ == '__main__':
    main()
