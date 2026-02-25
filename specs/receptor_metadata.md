# Receptor Metadata Structure

## Overview
Metadata files guide raw CIF → processed PDB transformation. Information is manually extracted from RCSB PDB database pages.

## Initial Metadata (Pre-Processing)

Created manually by reading PDB text files from RCSB database.

### Format
```json
{
    "target_ligand_resn": "FYP",
    "target_ligand_chain": "Z",
    "chains": {
        "A": "b2",
        "B": "a1",
        "C": "b2",
        "D": "a1_bzd",
        "E": "y2_bzd"
    }
}
```

### Fields

**`target_ligand_resn`** (string | null)
- Ligand residue name from PDB (e.g., "FYP", "DZP")
- `null` for apo structures

**`target_ligand_chain`** (string)
- Output chain for ligand after processing
- Always `"Z"` by convention

**`chains`** (object)
- Maps PDB chain IDs → subunit descriptions
- Only listed chains are KEPT (antibodies, extra chains deleted)
- Subunit naming:
  - `"b2"` = Beta-2
  - `"a1"` = Alpha-1
  - `"a2"` = Alpha-2
  - `"y2"` = Gamma-2
  - Append `"_bzd"` for chains at BZD binding interface (alpha + gamma)

## Workflow

### Step 1: receptor_prep.py
**Input:** Raw CIF + initial metadata
**Output:** `{pdb}_ligand.pdb` or `{pdb}_apo.pdb`

Process:
1. Load CIF
2. Keep only chains from `chains` dict (deletes antibodies F-I, water, NAG sugars, etc.)
3. If ligand present:
   - Find all copies of `target_ligand_resn`
   - Calculate BZD pocket center using `_bzd` chains
   - Select ligand closest to pocket
   - Move to chain Z
   - Delete other ligand copies
4. Save PDB

**NAG** = N-acetylglucosamine (sugar modification, automatically removed)

### Step 2: align_model_space.py (Master)
**Input:** `6X3U_ligand.pdb` + metadata
**Output:** `6X3U_aligned.pdb` + updated metadata

Process:
1. Center native ligand (FYP) at (0,0,0)
2. Standardize geometry (pore→Z, receptor COM→-X)
3. Ligand stays at chain Z
4. No copying (this is the master)

**Metadata updated:**
```json
{
    "chains": {...},
    "receptor_center_offset": [-26.978, 0.00003, -19.656],
    "alignment_reference": "6X3U",
    "ligands": {
        "Z": {
            "resname": "FYP",
            "source": "native",
            "native": true,
            "position": [0.0, 0.0, 0.0]
        }
    }
}
```

### Step 3: align_model_space.py (Holo→Reference)
**Input:** `6X3X_ligand.pdb` + `6X3U_aligned.pdb` (reference)
**Output:** `6X3X_aligned.pdb`

Process:
1. Align 6X3X receptor to 6X3U reference frame
2. Native DZP transforms with receptor
3. **Parse reference metadata** → read `ligands` dict
4. **Copy ligands from reference**:
   - Copy chain Z from 6X3U (FYP)
   - Bump native DZP from Z → Y
5. Save

**Result:** Z=FYP (copied), Y=DZP (native, bumped)

### Step 4: align_model_space.py (Apo→Reference)
**Input:** `9CRS_apo.pdb` + `6X3X_aligned.pdb` (reference with BOTH ligands)
**Output:** `9CRS_aligned.pdb`

Process:
1. Align 9CRS receptor to reference frame
2. **Parse reference metadata** → read `ligands` dict (has Z and Y)
3. **Copy ligands from reference**:
   - Copy chain Z (FYP)
   - Copy chain Y (DZP)
4. Save

**Result:** Z=FYP, Y=DZP (both copied)

## Key Insight

**Ligand chains grow cumulatively:**
- 6X3U master: Z only (FYP native)
- 6X3X aligned: Z, Y (FYP copied, DZP native bumped)
- 9CRS aligned: Z, Y (both copied)

**The copying logic:**
```python
# Read reference metadata
ref_meta = json.load("6X3X_aligned.json")

# Copy each ligand chain listed
for chain_id, lig_info in ref_meta["ligands"].items():
    if not lig_info["native"]:  # Don't re-copy native ligands
        copy_chain(ref_pdb, chain_id, target_pdb)
```

## Metadata Evolution

**Initial (manual):**
```json
{"target_ligand_resn": "FYP", "target_ligand_chain": "Z", "chains": {...}}
```

**After alignment:**
```json
{
    "chains": {...},
    "receptor_center_offset": [...],
    "alignment_reference": "6X3U",
    "ligands": {
        "Z": {"resname": "FYP", "source": "native", "native": true, "position": [0,0,0]},
        "Y": {"resname": "DZP", "source": "6X3X", "native": false, "position": [0.1,-0.2,0.05]}
    }
}
```

**Keys removed:** `target_ligand_resn`, `target_ligand_chain` (replaced by `ligands` dict)

## Chain Assignment Rules

- **A-E**: Receptor protein chains (pentamer)
- **Z**: First/reference ligand (FYP from master)
- **Y**: Second ligand (DZP, native or copied)
- **X, W, V...**: Additional ligands if needed (backwards from Z)
- **F-I**: Antibodies (deleted during prep)

## Data Sources

All chain/ligand info from RCSB PDB text files:
- `6X3U.txt`: "Entity ID" sections list chains + subunits
- "Ligands" section lists resnames (FYP, DZP, NAG, etc.)
- Manually identify receptor chains vs antibodies
- Manually map chains to subunit types from entity descriptions
