# Native Pose Validation Protocol

## Overview
When a crystal structure with a bound ligand is available, we can validate docking accuracy by calculating RMSD (Root Mean Square Deviation) between the docked poses and the native crystal structure pose. This is the gold standard for benchmarking molecular docking performance.

## Background
**Self-docking** (also called **re-docking**) is the process of:
1. Taking a ligand from its crystal structure
2. Removing it from the binding site
3. Docking it back into the same receptor
4. Comparing docked poses to the original crystal pose

This validates whether the docking algorithm can reproduce the experimentally observed binding mode.

## Success Criteria
Standard RMSD cutoffs for docking validation:
- **< 1.0 Å**: EXCELLENT - Highly accurate reproduction of native pose
- **< 2.0 Å**: GOOD - Acceptable binding mode prediction
- **> 2.0 Å**: POOR - Failed to reproduce native binding mode

## Our Workflow

### 1. Native Ligand Reference
For 6X3U (GABA_A α1β3γ2 with flumazenil):
- **Native ligand**: FYP (flumazenil) in chain Z
- **Source file**: `data/receptors/6X3U_ligand.pdb` (extracted from crystal structure)
- **Reference conformer**: `data/ligands/flumazenil.sdf` (RDKit-generated with proper bond info)

**Why we use the SDF instead of PDB:**
- Crystal structure PDBs lack explicit bond information
- RDKit cannot infer connectivity from coordinates alone for small molecules
- The prepared SDF has the same 3D coordinates but includes bond topology
- This is standard practice in computational chemistry

### 2. Docked Poses
GNINA outputs multiple poses ranked by score in a single SDF file:
- **File**: `data/output/6X3U_apo_flumazenil_out.sdf`
- **Format**: Multi-conformer SDF (9 poses)
- Each pose has properties: `minimizedAffinity`, `CNNscore`, `CNNaffinity`, `CNNvariance`

### 3. RMSD Calculation
**Script**: `scripts/pipeline/score_native_rmsd.py`

**Method**:
```python
from rdkit.Chem import rdMolAlign
rmsd = rdMolAlign.GetBestRMS(docked_mol, native_mol)
```

**Key features**:
- Uses `GetBestRMS()` which handles molecular symmetry automatically
- Aligns heavy atoms only (ignores hydrogens)
- Finds optimal superposition to minimize RMSD
- Robust to atom ordering differences

**Usage**:
```bash
conda activate docking
python scripts/pipeline/score_native_rmsd.py \
    --native data/ligands/flumazenil.sdf \
    --docked data/output/6X3U_apo_flumazenil_out.sdf \
    --out data/output/flumazenil_rmsd.txt
```

### 4. Results Interpretation
The script outputs:
1. **Per-pose RMSD**: RMSD for each docked conformer vs native
2. **Quality assessment**: EXCELLENT/GOOD/POOR classification
3. **Best pose**: Lowest RMSD across all poses
4. **Atom alignment**: Number of heavy atoms used in calculation

## 6X3U Flumazenil Results

### Summary
All 9 poses achieved **EXCELLENT** accuracy (< 1.0 Å RMSD):

| Pose | RMSD (Å) | Quality   |
|------|----------|-----------|
| 1    | 0.611    | EXCELLENT |
| 2    | 0.518    | EXCELLENT |
| 3    | 0.611    | EXCELLENT |
| **4**    | **0.478**    | **EXCELLENT** |
| 5    | 0.481    | EXCELLENT |
| 6    | 0.611    | EXCELLENT |
| 7    | 0.506    | EXCELLENT |
| 8    | 0.483    | EXCELLENT |
| 9    | 0.611    | EXCELLENT |

**Best pose: #4 with 0.478 Å RMSD**

### Interpretation
- **0.478 Å RMSD** is exceptionally accurate (typical benchmark is < 2.0 Å)
- All 9 top-ranked poses reproduced the native binding mode
- This validates our docking pipeline is working correctly
- The protocol (receptor prep, box definition, GNINA parameters) is sound
- We can confidently apply this to novel ligands (diazepam, docosanyl ferulate, etc.)

## Standard Practice
This validation approach follows established protocols in the literature:
- **Self-docking** is routine for validating docking setups
- RMSD to native pose is the standard metric (used in benchmarks like DUD, CASF, etc.)
- RDKit's `GetBestRMS()` handles symmetry (e.g., phenyl ring flips)
- Multi-pose analysis is critical (not just top-ranked pose)

## Files
- **Script**: `scripts/pipeline/score_native_rmsd.py`
- **Results**: `data/output/flumazenil_rmsd.txt`
- **Native reference**: `data/ligands/flumazenil.sdf`
- **Docked poses**: `data/output/6X3U_apo_flumazenil_out.sdf`

## Next Steps
1. Apply same validation to diazepam (if native pose available in another structure)
2. Use this confidence to dock novel compounds
3. For compounds without native poses, use other metrics:
   - Docking scores (minimizedAffinity, CNNscore)
   - Key residue interactions (His102, Phe77 distances)
   - Protein-ligand interaction fingerprints
   - Visual inspection in PyMOL

## References
- Bursulaya et al. (2003) - "Comparative study of several algorithms for flexible ligand docking"
- Mysinger et al. (2012) - DUD-E benchmark dataset
- Li et al. (2018) - "Comparative assessment of scoring functions"
- GNINA documentation on pose prediction accuracy
