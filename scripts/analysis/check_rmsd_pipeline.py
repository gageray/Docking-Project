#!/usr/bin/env python3
"""
Calculate RMSD between original aligned PDB and final PDBQT.
This verifies if coordinates are truly preserved through the pipeline.
"""

import numpy as np
from Bio.PDB import PDBParser, Superimposer

def get_ca_coords(pdb_path, chain_filter=None):
    """Extract CA atom coordinates from PDB/PDBQT file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)

    coords = []
    residue_ids = []

    for model in structure:
        for chain in model:
            if chain_filter and chain.id not in chain_filter:
                continue
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].coord)
                    residue_ids.append(f"{chain.id}:{residue.resname}{residue.id[1]}")

    return np.array(coords), residue_ids

print("="*60)
print("RMSD CALCULATION: PIPELINE COORDINATE PRESERVATION")
print("="*60)

# BZD chains only (D and E)
bzd_chains = ['D', 'E']

# Step 1: Original aligned PDB (BZD chains)
print("\n1. Loading original aligned PDB (chains D+E)...")
coords_original, res_ids_original = get_ca_coords("data/receptors/6X3U_aligned.pdb", chain_filter=bzd_chains)
print(f"   Found {len(coords_original)} CA atoms")

# Step 2: PyMOL stripped output
print("\n2. Loading PyMOL stripped PDB...")
coords_stripped, res_ids_stripped = get_ca_coords("data/receptors/prepped/6X3U_aligned_stripped_raw.pdb")
print(f"   Found {len(coords_stripped)} CA atoms")

# Step 3: PDB2PQR output
print("\n3. Loading PDB2PQR output...")
coords_pdb2pqr, res_ids_pdb2pqr = get_ca_coords("data/receptors/prepped/6X3U_aligned_stripped.pdb")
print(f"   Found {len(coords_pdb2pqr)} CA atoms")

# Step 4: Meeko PDBQT output
print("\n4. Loading Meeko PDBQT output...")
coords_pdbqt, res_ids_pdbqt = get_ca_coords("data/receptors/prepped/6X3U_aligned_apo.pdbqt")
print(f"   Found {len(coords_pdbqt)} CA atoms")

# Calculate RMSDs (no alignment, just raw difference)
print(f"\n{'='*60}")
print("RAW RMSD (NO SUPERPOSITION)")
print(f"{'='*60}")

if len(coords_original) == len(coords_stripped):
    diff = coords_original - coords_stripped
    rmsd_strip = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    print(f"\nOriginal → PyMOL Strip: {rmsd_strip:.6f} Å")
else:
    print(f"\nOriginal → PyMOL Strip: ATOM COUNT MISMATCH ({len(coords_original)} vs {len(coords_stripped)})")

if len(coords_stripped) == len(coords_pdb2pqr):
    diff = coords_stripped - coords_pdb2pqr
    rmsd_pdb2pqr = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    print(f"PyMOL Strip → PDB2PQR: {rmsd_pdb2pqr:.6f} Å")
else:
    print(f"PyMOL Strip → PDB2PQR: ATOM COUNT MISMATCH")

if len(coords_pdb2pqr) == len(coords_pdbqt):
    diff = coords_pdb2pqr - coords_pdbqt
    rmsd_meeko = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    print(f"PDB2PQR → Meeko PDBQT: {rmsd_meeko:.6f} Å")
else:
    print(f"PDB2PQR → Meeko PDBQT: ATOM COUNT MISMATCH")

if len(coords_original) == len(coords_pdbqt):
    diff = coords_original - coords_pdbqt
    rmsd_total = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    print(f"\n{'='*60}")
    print(f"TOTAL: Original → Final PDBQT: {rmsd_total:.6f} Å")
    print(f"{'='*60}")

    if rmsd_total < 0.001:
        print("\n✓ SUCCESS: Coordinates perfectly preserved (RMSD < 0.001 Å)")
    elif rmsd_total < 0.5:
        print(f"\n✓ ACCEPTABLE: Minor rounding errors only (RMSD = {rmsd_total:.3f} Å)")
    else:
        print(f"\n✗ FAILURE: Significant coordinate shift detected (RMSD = {rmsd_total:.3f} Å)")
else:
    print(f"\nOriginal → Final PDBQT: ATOM COUNT MISMATCH ({len(coords_original)} vs {len(coords_pdbqt)})")

print()
