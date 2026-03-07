# QM Ligand Preparation Script Documentation

## Location
`scripts/ligand/qm_ligand_prep.py`

---

## Purpose
Generate QM-optimized 3D conformers for molecular docking using GFN2-xTB via CREST. Integrates GNN-based pKa prediction to automatically determine correct protonation state at target pH.

---

## Workflow

### 1. Protonation State Determination (GNN)
- Takes input SMILES (any protonation state)
- Uses Moitessier GNN pKa predictor to analyze molecular context
- Returns dominant microstate SMILES at target pH (default 7.4)
- Calculates formal charge for CREST `--chrg` flag
- Falls back to input SMILES if GNN unavailable

### 2. 3D Seed Generation
- Converts corrected SMILES to 3D using RDKit ETKDG
- MMFF pre-optimization for reasonable starting geometry
- Exports to XYZ format for CREST

### 3. CREST Conformational Search
- Runs GFN2-xTB iMTD-GC conformer search
- Passes formal charge from GNN prediction
- Uses `--quick` mode for speed (minutes vs hours)
- Energy window: 6.0 kcal/mol (CREST default)
- Parses `crest_best.xyz` output

### 4. SDF Export
- Converts best conformer back to SDF
- Preserves topology from corrected SMILES
- Outputs to `data/ligands/{name}_qm.sdf`

---

## Modes

### Standard Mode
- Full conformational search on entire molecule
- No geometric constraints
- Use for rigid/small molecules

### Constrained Mode
- Identifies longest aliphatic chain via DFS
- Constrains chain to straight (all-trans) geometry
- Useful for long-chain molecules (docosanyl ferulate)
- Generates CREST `.cinp` constraint file
- Prevents chain folding during optimization

---

## CLI Usage

### Basic Example
```bash
python scripts/ligand/qm_ligand_prep.py \
  -n flumazenil \
  -s "CCOC(=O)C1=C2CN(C(=O)C3=C(N2C=N1)C=CC(=C3)F)C" \
  -t 4 \
  -m standard
```

### With pH Control
```bash
python scripts/ligand/qm_ligand_prep.py \
  -n ferulic_acid \
  -s "COC1=C(C=CC(=C1)/C=C/C(=O)O)O" \
  --ph 7.4 \
  -t 4
```

### Constrained Mode (Long Chains)
```bash
python scripts/ligand/qm_ligand_prep.py \
  -n docosanyl_ferulate \
  -s "CCCCCCCCCCCCCCCCCCCCCCOC(=O)/C=C/C1=CC(=C(C=C1)O)OC" \
  -m constrained \
  -t 4
```

### From Config File
```bash
python scripts/ligand/qm_ligand_prep.py -c config/flumazenil_job.json
```

**Config format:**
```json
{
  "name": "flumazenil",
  "smiles": "CCOC(=O)C1=C2CN(C(=O)C3=C(N2C=N1)C=CC(=C3)F)C",
  "outdir": "data/ligands",
  "threads": 4,
  "mode": "standard",
  "ph": 7.4
}
```

---

## Arguments

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-n, --name` | str | required | Ligand name |
| `-s, --smiles` | str | required | Input SMILES (neutral preferred) |
| `-o, --outdir` | str | `data/ligands` | Output directory |
| `-t, --threads` | int | 4 | CPU threads for CREST |
| `-m, --mode` | str | `standard` | `standard` or `constrained` |
| `--ph` | float | 7.4 | Target pH for protonation |
| `-c, --config` | str | None | JSON config (overrides CLI) |

---

## Output

### Files Generated

**Standard Mode:**
```
data/ligands/
├── flumazenil_crest/         # Working directory
│   ├── seed.xyz              # Input geometry
│   ├── crest_best.xyz        # Best conformer
│   └── crest_conformers.xyz  # All conformers
└── flumazenil_qm.sdf         # Final output
```

**Constrained Mode:**
```
data/ligands/
├── docosanyl_ferulate_crest/
│   ├── seed.xyz
│   ├── constrain.cinp        # CREST constraint file
│   ├── crest_best.xyz
│   └── crest_conformers.xyz
└── docosanyl_ferulate_qm.sdf
```

### SDF Properties
- 3D coordinates from GFN2-xTB optimization
- Correct protonation state at target pH
- Formal charge embedded in structure
- Heavy atoms + explicit hydrogens

---

## CREST Parameters

### Default Command
```bash
crest seed.xyz \
  --gfn2 \             # GFN2-xTB Hamiltonian
  --chrg +1 \          # Formal charge (from GNN)
  --ewin 6.0 \         # Energy window (kcal/mol)
  --quick \            # Fast search mode
  -T 4                 # CPU threads
```

### Search Mode: `--quick`
- Shortened MTD trajectory lengths
- Relaxed optimization convergence
- Capped iteration cycles
- **Runtime**: 5-15 minutes for drug-like molecules
- **vs Default**: Hours for exhaustive search

### Energy Window: 6.0 kcal/mol
- CREST default (NOT 40.0 like initial validation scripts)
- Retains conformers within thermal energy range
- Prevents exponential explosion of high-energy junk
- **DO NOT** increase to 40.0 for production

---

## GNN pKa Prediction Integration

### Library
`pka_predictor_moitessier` - Graph neural network trained on experimental pKa data

### Function: `get_correct_protonation_state()`

**Inputs:**
- `input_smiles`: Any protonation state
- `target_ph`: Target pH (default 7.4)
- `name`: Ligand name (logging)

**Outputs:**
- SMILES of dominant microstate at pH
- Formal charge for CREST

**Behavior:**
```python
# Example: Flumazenil
input_smiles = "CCOC(=O)C1=C2CN(C(=O)C3=C(N2C=N1)C=CC(=C3)F)C"  # Neutral
target_ph = 7.4

# GNN determines imidazole nitrogen should be protonated
corrected_smiles = "CCOC(=O)c1ncn2c1C[N@@H+](C)C(=O)c1cc(F)ccc1-2"  # Charged
formal_charge = +1  # Passed to CREST --chrg flag
```

**Fallback:**
If GNN unavailable or fails:
- Uses input SMILES as-is
- Calculates formal charge from input
- Prints warning

---

## Why This Approach

### Problem with Previous Pipeline
1. **Dimorphite-DL**: Rule-based protonation (pKa lookup tables)
2. **No context**: Treats all nitrogens the same
3. **Flumazenil failure**: Neutral imidazole incorrectly left neutral

### GNN Solution
1. **Whole-molecule context**: Analyzes electronic environment
2. **Conjugation aware**: Understands aromatic systems
3. **Flumazenil correct**: Imidazole protonated at pH 7.4

### MMFF → GFN2-xTB
1. **MMFF failure**: 99/100 conformers broke for charged species
2. **GFN2-xTB**: Quantum mechanical, handles charges natively
3. **Trade-off**: Minutes instead of milliseconds (acceptable)

---

## Validation Results

### Flumazenil (Before GNN + CREST)
- **Method**: Dimorphite-DL + ETKDG + MMFF
- **Valid conformers**: 1/100
- **RMSD to crystal**: 2.47 Å (failed threshold)
- **Issue**: Wrong protonation + force field failure

### Flumazenil (After GNN + CREST)
- **Method**: GNN + CREST/GFN2-xTB
- **Valid conformers**: All generated conformers valid
- **RMSD to crystal**: TBD (needs validation run)
- **Expected**: <2.0 Å (recovers crystal pose)

---

## Integration with Docking Pipeline

### Workflow
```
1. Input: Neutral SMILES from config.json
2. GNN pKa prediction → Charged SMILES at pH 7.4
3. CREST conformer generation → QM-optimized 3D geometry
4. Output: {name}_qm.sdf
5. GNINA docking with correct geometry + charge
```

### Config Updates Needed
Replace `ligand_prep` section in `config.json`:

**Before:**
```json
"ligand_prep": {
  "output_dir": "data/ligands",
  "ligands": {
    "flumazenil": {
      "smiles": "CCOC(=O)C1=C2CN(C(=O)C3=C(N2C=N1)C=CC(=C3)F)C",
      "mode": "single"
    }
  }
}
```

**After:**
```json
"qm_ligand_prep": {
  "output_dir": "data/ligands",
  "threads": 4,
  "ph": 7.4,
  "ligands": {
    "flumazenil": {
      "smiles": "CCOC(=O)C1=C2CN(C(=O)C3=C(N2C=N1)C=CC(=C3)F)C",
      "mode": "standard"
    },
    "docosanyl_ferulate": {
      "smiles": "CCCCCCCCCCCCCCCCCCCCCCOC(=O)/C=C/C1=CC(=C(C=C1)O)OC",
      "mode": "constrained"
    }
  }
}
```

---

## Dependencies

- **RDKit**: SMILES parsing, ETKDG, MMFF
- **CREST**: GFN2-xTB conformer search (conda-forge)
- **pka_predictor_moitessier**: GNN pKa prediction (pip)
- **Python 3.8+**

---

## Common Issues

### GNN Import Error
**Error**: `ImportError: No module named 'pka_predictor_moitessier'`
**Fix**: Install GNN library or script falls back to input SMILES

### CREST Not Found
**Error**: `CREST binary not found in PATH`
**Fix**: `conda install -c conda-forge crest`

### ETKDG Embedding Fails
**Behavior**: Script tries bare minimum embedding as fallback
**Rare**: Usually only fails for complex macrocycles

### Constrained Mode No Tail Found
**Behavior**: Falls back to standard mode automatically
**When**: Molecule has no aliphatic chain >4 carbons

---

## Performance

| Molecule Size | Mode | Threads | Runtime |
|---------------|------|---------|---------|
| 20-30 atoms (flumazenil) | standard | 4 | 5-10 min |
| 30-50 atoms (ferulate) | standard | 4 | 10-20 min |
| 50+ atoms (DF) | constrained | 4 | 15-30 min |

**Bottleneck**: CREST conformational search (not GNN prediction)

---

## Future Enhancements

1. **Batch mode**: Process multiple ligands in parallel
2. **Ensemble output**: Save all CREST conformers (not just best)
3. **Energy filtering**: Custom `--ewin` per ligand
4. **Validation integration**: Auto-compare to crystal structures
