#!/usr/bin/env python3
"""
Check if center of mass shifts during receptor preparation pipeline.
"""

import numpy as np
from Bio.PDB import PDBParser

def calculate_com(pdb_path, chain_filter=None):
    """Calculate center of mass for a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)

    coords = []
    for model in structure:
        for chain in model:
            if chain_filter and chain.id not in chain_filter:
                continue
            for residue in chain:
                for atom in residue:
                    coords.append(atom.coord)

    if not coords:
        return None

    coords = np.array(coords)
    com = np.mean(coords, axis=0)
    return com

# Check COM at each stage
print("="*60)
print("CENTER OF MASS TRACKING")
print("="*60)

# Original aligned PDB (all chains)
com_original_all = calculate_com("data/receptors/6X3U_aligned.pdb")
print(f"\n1. Original aligned PDB (all chains):")
print(f"   COM: {com_original_all}")

# Original aligned PDB (only BZD chains D and E)
com_original_bzd = calculate_com("data/receptors/6X3U_aligned.pdb", chain_filter=['D', 'E'])
print(f"\n2. Original aligned PDB (BZD chains D+E only):")
print(f"   COM: {com_original_bzd}")

# PyMOL stripped output (before PDB2PQR)
com_stripped = calculate_com("data/receptors/prepped/6X3U_aligned_stripped_raw.pdb")
print(f"\n3. After PyMOL strip (raw):")
print(f"   COM: {com_stripped}")

# After PDB2PQR protonation
com_protonated = calculate_com("data/receptors/prepped/6X3U_aligned_stripped.pdb")
print(f"\n4. After PDB2PQR protonation:")
print(f"   COM: {com_protonated}")

# Calculate shifts
if com_original_bzd is not None and com_stripped is not None:
    shift_pymol = com_stripped - com_original_bzd
    shift_magnitude = np.linalg.norm(shift_pymol)

    print(f"\n{'='*60}")
    print("SHIFT ANALYSIS")
    print(f"{'='*60}")
    print(f"PyMOL strip shift: {shift_pymol}")
    print(f"Shift magnitude: {shift_magnitude:.6f} Å")

    if shift_magnitude < 0.001:
        print("✓ NO SHIFT - Coordinates preserved")
    else:
        print(f"✗ SHIFT DETECTED - Coordinates moved {shift_magnitude:.3f} Å")

print(f"{'='*60}\n")
