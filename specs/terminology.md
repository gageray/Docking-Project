# Terminology & Conventions

## Receptor States

**Holo Structure**
- Receptor with ligand bound in crystal structure
- Examples: 6X3U (with FYP), 6X3X (with DZP)
- Used as alignment references
- Ligand position validates pocket geometry

**Apo Structure**
- Receptor without ligand (empty binding site)
- Examples: 9CRS, 9CSB
- Must be aligned to holo reference structure
- Validation ligands copied from reference

## Ligand Types

**FYP** - Flumazenil
- Benzodiazepine antagonist (blocks binding)
- Reference ligand from 6X3U
- Always at chain Z in aligned structures

**DZP** - Diazepam
- Benzodiazepine agonist (activates receptor)
- Validation ligand from 6X3X
- Typically at chain Y in aligned structures

**FA** - Ferulic Acid
- Natural compound, novel docking target
- Small, rigid structure

**DF** - Docosanyl Ferulate
- FA with 22-carbon tail
- Requires special conformer handling

## Receptor Subtypes

**α1β2γ2** (Alpha-1)
- 6X3U, 6X3X, 9CRS
- Most common GABA_A subtype
- High benzodiazepine affinity

**α2β2β3γ2** (Alpha-2)
- 9CSB
- Different pharmacology than α1
- Used for selectivity comparison

## Chain Naming

**Protein Chains**
- A-E: Receptor subunits (pentamer)
- F-I: Antibodies (removed during prep)

**Ligand Chains**
- Z: Primary/reference ligand (FYP in all aligned)
- Y: Secondary ligand (DZP, or native bumped from Z)
- X, W, V...: Additional ligands (backwards from Z)

## Post-Translational Modifications (PTMs)

**NAG** - N-Acetylglucosamine
- Sugar modification on protein surface
- Adds conformational flexibility
- Standard practice: **REMOVE for docking**
  - Not in binding pocket
  - Increases complexity without benefit
  - Docking software often strips automatically
- Keep for MD simulations (affects dynamics)

**Other PTMs**
- Phosphorylation, methylation, etc.
- Usually absent in cryo-EM structures
- If present, evaluate case-by-case

## Processing Stages

**Raw CIF**
- Direct from PDB database
- Contains: receptor, antibodies, ligands, water, sugars
- Example: `6X3U.cif`

**Prepped PDB**
- After receptor_prep.py
- Kept: receptor chains (A-E) + target ligand (→Z)
- Removed: antibodies, extra ligands, water, NAG
- Example: `6X3U_ligand.pdb` or `9CRS_apo.pdb`

**Aligned PDB**
- After align_model_space.py
- Pocket centered at (0, 0, 0)
- Pore along Z-axis
- Validation ligands added
- Example: `6X3U_aligned.pdb`

## Coordinate System

**Origin** (0, 0, 0)
- Center of BZD binding pocket
- Reference point for all measurements

**Z-axis**
- Ion channel pore direction
- Points "up" through membrane

**-X axis**
- Receptor center of mass
- Pocket faces +X direction

## Alignment Strategy

**Master Reference** - 6X3U
- First structure aligned (independent)
- Defines coordinate frame for all others
- FYP centered at origin

**Sequential Alignment**
- 6X3U (master) → 6X3X (holo) → 9CRS/9CSB (apo)
- Each inherits ligands from previous
- Cumulative ligand collection: Z, then Z+Y

## Metadata Evolution

**Initial** (pre-processing)
```json
{
  "target_ligand_resn": "FYP",
  "target_ligand_chain": "Z",
  "chains": {"A": "b2", ...}
}
```

**After Alignment**
```json
{
  "chains": {...},
  "receptor_center_offset": [...],
  "alignment_reference": "6X3U",
  "ligands": {
    "Z": {"resname": "FYP", "native": true, "position": [0,0,0]},
    "Y": {"resname": "DZP", "native": false, "position": [...]}
  }
}
```

## Docking Terminology

**Binding Affinity**
- GNINA score in kcal/mol
- More negative = stronger binding
- Typical range: -6 to -12 kcal/mol

**RMSD** (Root Mean Square Deviation)
- Measure of pose similarity in Ångströms
- For re-docking validation: RMSD < 2.0 Å = success
- For cross-comparison: lower RMSD = more similar binding mode

**Pocket**
- Benzodiazepine binding site
- Interface between α (or α2) and γ2 subunits
- ~10 Å radius from origin

**BZD Chains**
- Chains forming the binding pocket
- Tagged with `_bzd` in metadata
- Used for alignment in apo mode

## File Naming Conventions

**PDB Downloads**
- `{PDB_ID}.cif` - Raw structure
- `{PDB_ID}.txt` - RCSB metadata

**Processed Structures**
- `{PDB_ID}_ligand.pdb` - Holo, after prep
- `{PDB_ID}_apo.pdb` - Apo, after prep
- `{PDB_ID}_aligned.pdb` - After alignment

**Metadata**
- `{pdb_id}_metadata.json` - Initial/input
- `{pdb_id}_aligned.json` - After alignment

**Validation Ligands**
- Extracted ligands stored in `validation_ligands/`
- Format: `{RESNAME}_{fullname}.pdb`
- Only needed if extracting separately (currently not used)

## Standard Practices

**Why Remove Water?**
- Too many molecules (thousands)
- Adds no structural information for rigid docking
- Solvent effects handled implicitly

**Why Remove Antibodies?**
- Used for crystallization/cryo-EM, not physiological
- Block access to regions of interest
- Add unnecessary atoms

**Why Keep PTMs in MD but not Docking?**
- **MD**: Long timescale, flexibility matters
- **Docking**: Static pose, PTMs add noise not signal
