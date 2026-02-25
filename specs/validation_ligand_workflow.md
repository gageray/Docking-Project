# Validation Ligand Injection & Metrics

## Overview
All receptors aligned to common reference frame (6X3U) with validation ligands (FYP, DZP) injected at standardized positions. Enables visual and computational comparison of docking results across receptor subtypes.

## Validation Strategy
1. **Re-dock known ligands** (FYP, DZP) into their crystal structures → measure RMSD to validate pipeline
2. **Cross-dock** FYP and DZP into all receptors → compare binding modes
3. **Dock novel ligands** (FA, DF) → compare interactions to benzodiazepine reference poses

---

## Workflow

### **Master Reference: 6X3U (FYP)**

### **Step 1: Align Master Receptor**
```bash
align_model_space.py \
  --pdb 6X3U_ligand.pdb \
  --json 6x3u_metadata.json \
  --out_pdb 6X3U_aligned.pdb \
  --out_json 6x3u_aligned_metadata.json \
  --master_receptor 6X3U
```
**Result:**
- Native FYP at chain Z at (0, 0, 0)
- Pore aligned to Z-axis
- Receptor COM on -X axis
- Independent alignment (no reference)

### **Step 2: Align Holo Receptor to Reference (6X3X)**
```bash
align_model_space.py \
  --pdb 6X3X_ligand.pdb \
  --json 6x3x_metadata.json \
  --out_pdb 6X3X_aligned.pdb \
  --out_json 6x3x_aligned_metadata.json \
  --ref_pdb 6X3U_aligned.pdb \
  --master_receptor 6X3U
```
**Result:**
- Receptor aligned to 6X3U reference frame
- Native DZP transformed with receptor
- Copy FYP from 6X3U chain Z → inject at chain Z
- Bump native DZP from Z → Y
- Final: **Z=FYP (copied), Y=DZP (native, transformed)**

### **Step 3: Align Apo Receptors (9CRS, 9CSB)**
```bash
align_model_space.py \
  --pdb 9CRS_apo.pdb \
  --json 9crs_metadata.json \
  --out_pdb 9CRS_aligned.pdb \
  --out_json 9crs_aligned_metadata.json \
  --ref_pdb 6X3U_aligned.pdb \
  --master_receptor 6X3U
```
**Result:**
- Receptor aligned to 6X3U reference frame
- Copy FYP from 6X3U chain Z → inject at chain Z
- Copy DZP from 6X3X chain Y → inject at chain Y
- Final: **Z=FYP (ref), Y=DZP (ref)**

---

## Chain Assignment Rules

### **Chain Naming Convention**
- Work backwards alphabetically from Z
- **Z**: Always reference ligand (FYP from 6X3U)
- **Y**: Second validation ligand (DZP from 6X3X) OR native ligand (if bumped)
- **X, W, V...**: Additional ligands as needed

### **Native Ligand Bumping**
When validation ligands injected into holo receptor:
1. Native ligand at Z → bumped to Y
2. Reference ligand injected at Z
3. Additional validation ligands at X, W, etc.

---

## Metadata Structure

```json
{
  "alignment_reference": "6X3U",
  "chains": {
    "A": "b2",
    "B": "a1",
    "C": "b2",
    "D": "a1_bzd",
    "E": "y2_bzd"
  },
  "receptor_center_offset": [-26.978, 0.0, -19.656],
  "ligands": {
    "Z": {
      "resname": "FYP",
      "source": "6X3U",
      "native": false,
      "position": [0.0, 0.0, 0.0]
    },
    "Y": {
      "resname": "DZP",
      "source": "native",
      "native": true,
      "position": [0.05, -0.02, 0.01]
    }
  },
  "alignment_quality": {
    "pocket_rmsd_vs_reference": 0.28,
    "pocket_radius": 10.0,
    "n_pocket_atoms": 342,
    "key_residues": {
      "D_His102": {
        "CA_position": [x, y, z],
        "distance_to_ligand_Z": 3.2
      },
      "E_Phe77": {
        "CA_position": [x, y, z],
        "distance_to_ligand_Z": 3.8
      }
    }
  }
}
```

---

## Validation Metrics

### **1. Alignment Quality (Receptor)**
**Pocket RMSD:**
- Select all heavy atoms within 10Å of pocket center (0,0,0)
- Superimpose vs 6X3U reference pocket
- RMSD < 0.5 Å = excellent alignment
- RMSD < 1.0 Å = acceptable alignment

**Key Residue Tracking:**
- CA positions of binding site residues (D_His102, E_Phe77, etc.)
- Distance between key residues
- Validates pocket geometry preservation

### **2. Ligand Positioning (Pre-Docking)**
**Ligand COM:**
- Record 3D position of each ligand centroid
- Reference ligands should have identical positions across all receptors
- Native ligands show positional variation due to COM differences

**Distance to Key Residues:**
- Min distance from ligand to D_His102
- Min distance from ligand to E_Phe77
- Baseline for comparing docked poses

### **3. Docking Validation (Post-Docking)**
**Ligand RMSD:**
- Re-dock FYP into 6X3U → RMSD vs crystal pose
- Re-dock DZP into 6X3X → RMSD vs crystal pose
- Success threshold: RMSD < 2.0 Å
- Good threshold: RMSD < 1.0 Å

**Binding Affinity:**
- GNINA docking scores (kcal/mol)
- Compare FYP vs DZP vs FA vs DF
- Compare across receptor subtypes (α1 vs α2)

**Interaction Analysis:**
- H-bond counts and geometry
- π-π stacking (aromatic residues)
- Distance to His102, Phe77, other key residues

---

## Implementation Details

### **align_model_space.py Logic**

**Case A: Master Receptor (no --ref_pdb)**
- Independent alignment
- Native ligand centered at origin
- No ligand injection

**Case B: Holo to Reference (has ligand AND --ref_pdb)**
- Load target receptor with native ligand
- Align to reference (superimpose pocket residues)
- Native ligand transforms with receptor
- Copy reference ligands from master
- Bump native ligand to next chain
- Inject reference ligands starting at Z

**Case C: Apo to Reference (no ligand AND --ref_pdb)**
- Align to reference
- Copy reference ligands from master
- Inject starting at Z

### **Pocket Selection**
- All heavy atoms within 10Å radius of (0, 0, 0)
- Excludes backbone atoms from non-BZD chains
- Focuses on binding site residues

### **RMSD Calculation**
```python
# Pocket RMSD
ref_pocket = select_pocket(ref_structure, radius=10.0)
target_pocket = select_pocket(target_structure, radius=10.0)
rmsd = superimpose_and_calc_rmsd(ref_pocket, target_pocket)

# Ligand RMSD (post-docking)
crystal_coords = get_ligand_coords(crystal_structure, "FYP")
docked_coords = get_ligand_coords(docked_pose, "FYP")
rmsd = calc_rmsd(crystal_coords, docked_coords)
```

---

## Validation Outputs

### **For Publication/Reporting**
1. **Alignment Quality Table:**
   - Receptor | Pocket RMSD | N atoms | Key residue distances

2. **Re-docking Validation:**
   - Ligand | Receptor | RMSD | Docking Score

3. **Cross-docking Results:**
   - Ligand | Receptor | Score | Key interactions

4. **Novel Ligand Analysis:**
   - FA/DF positioning relative to FYP/DZP
   - Interaction with His102/Phe77
   - Tail conformation (DF only)

### **Visualization**
- PyMOL overlays of validation ligands across receptors
- Docked poses overlaid with crystal structures
- Interaction maps (2D ligand-protein diagrams)

---

## File Naming Convention

**Aligned receptors:**
- `6X3U_aligned.pdb` - master reference
- `6X3X_aligned.pdb` - holo aligned to master
- `9CRS_aligned.pdb` - apo aligned to master

**Metadata:**
- `6x3u_aligned_metadata.json`
- `6x3x_aligned_metadata.json`
- `9crs_aligned_metadata.json`

**Validation ligands:**
- Extracted from aligned structures as needed
- Chain Z = FYP (reference)
- Chain Y = DZP (reference or native)

---

## Notes

- **Floating point precision:** Positions at "0,0,0" are actually ~1e-5 to 1e-8 due to floating point arithmetic
- **COM vs geometric center:** Using geometric center (unweighted mean) for small ligands; adequate precision
- **DF (docosanyl ferulate):** Always use "DF" abbreviation
- **Pocket radius:** 10Å standard for BZD binding site; captures all interacting residues
- **Reference frame:** All alignment relative to 6X3U (α1β2γ2 with FYP)
