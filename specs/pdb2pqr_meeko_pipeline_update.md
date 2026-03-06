# Receptor Preparation Pipeline Update (PDB2PQR + Meeko)

## Objective
Replaced the legacy OpenBabel-based receptor preparation workflow with a highly accurate, chemically robust pipeline utilizing **PDB2PQR** and **Meeko**.

## Key Workflow Changes

1. **Chain Stripping (PyMOL)**: 
   - Uses `cmd.remove("not polymer.protein")` during processing to strictly isolate the receptor protein and automatically eliminate cryo-EM artifacts (like extraneous FYP ligands or buffers).
   
2. **Protonation & Hydrogen Addition (PDB2PQR)**: 
   - Assigns chemically accurate protonation states at pH 7.4 using PropKa. 
   - Crucially, the `--nodebump` flag was added to lock heavy atoms precisely in their original coordinates. This ensures the structural alignment framework is perfectly preserved for GNINA without local residue rotation compensations.

3. **PDB Sanitization (`fix_pdb_elements`)**: 
   - A custom Python function (`fix_pdb_elements`) was implemented in the conversion script to restore standard PDB element symbols omitted by PDB2PQR in columns 77-78. 
   - This explicitly prevents downstream `RuntimeError: Post-condition Violation Element '' not found` failures when parsed by RDKit inside Meeko.

4. **PDBQT Conversion (Meeko)**: 
   - Uses `mk_prepare_receptor.py` natively using standard PyPI syntax flags (`-i`, `-o` without `.pdbqt` extension, and `-p`) to securely output the final `.pdbqt` file.

5. **Environment Configuration**: 
   - Permanently removed OpenBabel from the project requirements. 
   - Switched to utilizing `pdb2pqr` and `meeko` directly installed from the `conda-forge` channel. This has been updated in the worker's `setup_env.py` script.

## Verification Checkpoints Passed
- Confirmed there is zero coordinate shifting or rotation for any alpha or carbon backbones between the initial extraction and the final output.
- Validated that the final PDBQT output specifically isolates the correct BZD pocket chains (e.g. D and E) effectively matching the parameters listed in `metadata.json`.
