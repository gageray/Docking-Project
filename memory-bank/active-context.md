# Active Context

Currently implemented the automated data pipeline for docking.
- Configured JSON-based single source of truth (`config.json`) for the pipeline.
- `ligand_prep.py` constructed to run RDKit conformer generation, dimorphite-dl protonation, and deterministic grid searches for conformers.
- `receptor_prep.py` prepared to convert and clean CIF files to PDBs using OpenMM PDBFixer.
- `box_definition.py` created to run sequence-aware CA alignments across receptors and map docking bounding box centroids.
- All scripts are now housed in `scripts/`.
- Visualization Automation added via three PyMOL execution scripts:
  - `vis_ligands.py`: Validates distinct tail distributions vs headgroup alignment.
  - `vis_box.py`: Automatically drops PyMOL pseudoatoms onto computed target box centroids for overlap evaluation vs reference crystal.
  - `vis_results.py`: Post-docking template to analyze interactions and collisions between the target surface and the docked conformers.
- Environment setup and requirements defined in `requirements.txt` and `README.md` (now includes `pymol-open-source`).
- Refined pipeline scripts following code review: improved SMARTS matching, CA indexing, constraint logic, and added logging via a shared `utils.py`. PyMOL scripts parameterized and operator logic debugged.
