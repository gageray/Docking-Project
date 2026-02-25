# Validation Metrics Baseline Rationale

## Overview
The validation metrics script ([calc_validation_metrics.py](../scripts/receptor/calc_validation_metrics.py)) calculates pocket RMSD and ligand-to-residue distance measurements to validate receptor alignment quality and compare pocket geometry across receptor subtypes.

## Baseline Choice: Universal 6X3U FYP Reference

### Decision
All distance measurements (Δ values) are calculated relative to **6X3U FYP (chain Z)** as a universal baseline across all receptors.

### Rationale

#### Why NOT Native Ligand Baseline?
Initially considered using each receptor's native ligand as its own baseline:
- 6X3U: Compare validation ligands to native FYP
- 6X3X: Compare validation ligands to native DZP
- 9CRS/9CSB: No native (apo structures)

**Problems:**
1. **Inconsistent zero point**: What "Δ = 0" means changes per receptor
2. **Cross-receptor comparison impossible**: Can't compare α1 vs α2 geometry
3. **Apo structures have no native**: Would need special handling
4. **Confusing interpretation**: -1.6Å from FYP vs -1.6Å from DZP are different absolute positions

#### Why Universal 6X3U FYP Baseline?
**The alignment establishes a common reference frame:**
1. 6X3U is aligned first, placing FYP at origin (0,0,0)
2. All other receptors align their BZD pocket to 6X3U pocket geometry
3. Pocket RMSD confirms geometric alignment (~0.4Å)
4. **Distance measurements reflect pocket geometry, not alignment quality**

**Advantages:**
- **Consistent zero point**: Δ = 0 means "same distance as 6X3U FYP to His102/Phe77"
- **Cross-subtype comparison**: Direct comparison of α1 (His102) vs α2 (His101) geometry
- **Absolute meaning**: Positive Δ = further from residue, negative Δ = closer to residue
- **Works for apo**: No special handling needed for 9CRS/9CSB

### What the Metrics Measure

#### Pocket RMSD
- **Calculation**: CA atoms within 10Å of origin, no transformation applied
- **Meaning**: How well the binding pocket geometry aligns across receptors
- **Result**: ~0.4Å for all receptors → excellent alignment

#### Δ Ligand-Residue Distances
- **Baseline**: 6X3U FYP-His102 = 10.125Å, FYP-Phe77 = 7.282Å
- **Meaning**: How ligand position in this pocket compares to FYP position in reference pocket
- **NOT measuring**: Alignment quality (that's what RMSD is for)
- **IS measuring**: Absolute pocket geometry differences

### Key Insight: Pocket Geometry vs Alignment

**The pocket is well-aligned (RMSD ~0.4Å), but the residues are in different absolute positions.**

Example - 9CSB (α2 subtype):
- Pocket RMSD: 0.420Å (good alignment)
- Δ FYP-His: +3.025Å (His is 3Å further from FYP position than in 6X3U)

**Interpretation**: The BZD binding pocket aligns well structurally, but the **α2 subtype has His101 in a different position** than α1's His102. This geometric difference is **biologically meaningful** - it may explain differential FYP binding affinity between α1 and α2 subtypes.

### Results Interpretation

| Receptor | Subtype | Δ FYP-His | Interpretation |
|----------|---------|-----------|----------------|
| 6X3U | α1β2γ2 | +0.000Å | Reference baseline |
| 6X3X | α1β2γ2 | +0.158Å | Nearly identical α1 geometry |
| 9CRS | α1β2γ2 | -0.051Å | Nearly identical α1 geometry |
| 9CSB | α2β2β3γ2 | **+3.025Å** | **α2 has different His position** |

**Conclusion**: The +3Å difference in 9CSB shows that **α2 subtypes have different BZD pocket geometry** compared to α1, specifically in the position of the critical His residue at the α-γ interface.

## Implementation Details

### Residue Selection
- **His residue**: Uses reference His102 for ALL comparisons (even 9CSB which has His101)
  - **Why**: Maintains consistent baseline; measures to same pocket position
  - **Note**: 9CSB His101 is at a different absolute position than 6X3U His102
- **Phe residue**: Phe77 on γ2 chain (conserved across all receptors)

### Distance Calculation
- `cmd.distance(lig_sel, res_sel)` returns minimum distance between selections
- `lig_sel`: All heavy atoms of ligand (no hydrogens)
- `res_sel`: CA atom of target residue

### Output Format
```
Reference: 6X3U_aligned.pdb (FYP baseline at chain Z)
Baseline distances: FYP-His102=10.125Å, FYP-Phe77=7.282Å
---
Target | Pocket RMSD | CA Match | Native Lig | Δ FYP-His | Δ DZP-His | Δ FYP-Phe | Δ DZP-Phe
```

## Future Considerations

### For Docking Validation
After docking validation ligands (FYP/DZP) to each receptor:
- Compare docked pose RMSD to crystallographic position
- Calculate ligand-residue distances in docked poses
- Use same universal baseline for consistency

### For Additional Receptors
When adding new receptors:
- Align to 6X3U reference frame
- Use same FYP baseline for distance comparisons
- Note if receptor has different subtype (e.g., α3, α5) for interpretation

## See Also
- [Validation Ligand Workflow](validation_ligand_workflow.md) - Alignment and ligand copying strategy
- [Receptor Metadata](receptor_metadata.md) - Metadata structure and evolution
- [Terminology](terminology.md) - Definitions of holo/apo, subtypes, etc.
