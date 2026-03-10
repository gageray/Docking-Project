#!/usr/bin/env python3
"""
Check if center of mass shifts during receptor preparation pipeline.
Track each chain individually.
"""

import numpy as np
from Bio.PDB import PDBParser

def calculate_com_per_chain(pdb_path):
    """Calculate center of mass for each chain in a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)

    chain_coms = {}
    for model in structure:
        for chain in model:
            coords = []
            for residue in chain:
                for atom in residue:
                    coords.append(atom.coord)

            if coords:
                coords = np.array(coords)
                chain_coms[chain.id] = np.mean(coords, axis=0)

    return chain_coms

# Check COM at each stage
print("="*60)
print("CENTER OF MASS TRACKING (PER CHAIN)")
print("="*60)

# Original aligned PDB
print("\n1. ORIGINAL ALIGNED PDB (6X3U_aligned.pdb):")
com_original = calculate_com_per_chain("data/receptors/6X3U_aligned.pdb")
for chain_id, com in sorted(com_original.items()):
    print(f"   Chain {chain_id}: {com}")

# PyMOL stripped output (before PDB2PQR)
print("\n2. AFTER PYMOL STRIP (6X3U_aligned_stripped_raw.pdb):")
com_stripped = calculate_com_per_chain("data/receptors/prepped/6X3U_aligned_stripped_raw.pdb")
for chain_id, com in sorted(com_stripped.items()):
    print(f"   Chain {chain_id}: {com}")

# After PDB2PQR protonation
print("\n3. AFTER PDB2PQR (6X3U_aligned_stripped.pdb):")
com_protonated = calculate_com_per_chain("data/receptors/prepped/6X3U_aligned_stripped.pdb")
for chain_id, com in sorted(com_protonated.items()):
    print(f"   Chain {chain_id}: {com}")

# Calculate shifts for BZD chains
print(f"\n{'='*60}")
print("SHIFT ANALYSIS (BZD CHAINS D & E)")
print(f"{'='*60}")

for chain_id in ['D', 'E']:
    if chain_id in com_original and chain_id in com_stripped:
        shift = com_stripped[chain_id] - com_original[chain_id]
        magnitude = np.linalg.norm(shift)

        print(f"\nChain {chain_id}:")
        print(f"  Original COM: {com_original[chain_id]}")
        print(f"  Stripped COM: {com_stripped[chain_id]}")
        print(f"  Shift vector: {shift}")
        print(f"  Magnitude: {magnitude:.6f} Å")

        if magnitude < 0.001:
            print(f"  ✓ NO SHIFT - Coordinates preserved")
        else:
            print(f"  ✗ SHIFT DETECTED - Coordinates moved {magnitude:.3f} Å")

print(f"{'='*60}\n")
