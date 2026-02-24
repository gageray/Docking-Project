# Active Context

## Current State (Feb 23, 2026)

### Project Structure Reorganization
- Scripts reorganized into subdirectories:
  - `scripts/ligand/` - Ligand preparation pipeline
  - `scripts/receptor/` - Receptor preparation and alignment
  - `scripts/visualization/` - PyMOL rendering and validation
- Old flat scripts deleted from git (pending commit)
- New spec files added: `docking optimization.txt`, `vis_rotation.txt`

### Implemented Pipeline
1. **Ligand Prep** (`scripts/ligand/ligand_prep.py`)
   - RDKit conformer generation with Dimorphite-DL pH 7.4 protonation
   - Single conformer mode for simple ligands
   - Grid search mode for complex ligands (72 conformers)
   - Energy filtering (>20 kcal/mol cutoff)

2. **Receptor Prep** (`scripts/receptor/receptor_prep.py`)
   - CIF → PDB conversion via Biopython
   - Chain extraction and filtering
   - PDBFixer for missing atoms/residues/hydrogens
   - pH 7.4 protonation state

3. **Geometry Standardization** (`scripts/receptor/align_model_space.py`)
   - Aligns pore along Z-axis via PCA
   - Centers binding pocket at (0,0,0)
   - Orients receptor COM on -X axis
   - Handles holo mode (global PCA+translation) and apo mode (10A radial alignment against reference BZD chains)

4. **Box Definition** (`scripts/receptor/box_definition.py`)
   - Sequence-aware CA superimposition
   - Reference centroid mapping to targets
   - JSON export of box parameters

5. **Visualization Suite** (`scripts/visualization/`)
   - `vis_ligands.py` - Conformer distribution validation
   - `vis_box.py` - Binding box overlay visualization
   - `vis_results.py` - Post-docking analysis
   - `render_single.py` - Individual structure rendering
   - `render_sweep.py` - Batch rendering pipeline
   - `check_view.py` - Camera angle verification

### Supporting Infrastructure
- `scripts/utils.py` - Shared logging utilities
- Unit tests in `tests/test_ligand_prep.py`
- `.claude/claude.md` - Project documentation for Claude Code
- Updated `.gitignore` with comprehensive exclusions
- `.agent/` folder with rules and workflows
- `memory-bank/` for session context tracking

### Configuration
- `config.json` - Central pipeline configuration
- All paths, parameters, and targets defined in JSON
- No hardcoded constants in scripts

### Code Quality
- Post-review refinements applied
- Improved SMARTS matching for ester identification
- Fixed CA indexing in alignment
- Enhanced constraint logic for tail locking
- Comprehensive logging throughout
