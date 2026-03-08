#!/usr/bin/env python3
"""
QupKake-based pKa Prediction and Protonation State Correction
Uses QupKake QM calculations to determine protonation state at target pH.
Defaults: Tautomerization (-t) is enabled, multiprocessing (-mp) set to 4 cores.
"""

import subprocess
import tempfile
import os
from rdkit import Chem


def get_pka_and_state(input_smiles: str, target_ph: float = 7.4, name: str = "ligand"):
    """
    Executes QupKake to calculate whole-molecule QM/ML pKa, adjusts protonation
    states based on target_ph, and returns the corrected SMILES with pKa values.

    Args:
        input_smiles: Input SMILES (any protonation state)
        target_ph: Target pH (default 7.4)
        name: Ligand name for logging

    Returns:
        Tuple of (corrected_smiles, pka_records_dict)
        where pka_records_dict = {"Atom 5 (N)": 2.87, ...}
    """
    print(f"[{name}] Analyzing protonation state via QupKake at pH {target_ph}...")

    with tempfile.TemporaryDirectory() as tmpdir:
        # QupKake uses --root parameter and saves to {root}/output/qupkake_output.sdf
        out_sdf = os.path.join(tmpdir, "output", "qupkake_output.sdf")

        try:
            # Execute QupKake with --root pointing to tmpdir
            # Added -t (tautomerize) and -mp 4 (multiprocessing with 4 cores) by default
            subprocess.run(
                ["qupkake", "smiles", input_smiles, "--root", tmpdir, "-n", name, "-t", "-mp", "4"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
        except subprocess.CalledProcessError as e:
            print(f"[{name}] ERROR: QupKake execution failed: {e}")
            print(f"[{name}] Falling back to input SMILES")
            return input_smiles, {}
        except FileNotFoundError:
            print(f"[{name}] ERROR: QupKake not found in PATH")
            print(f"[{name}] Falling back to input SMILES")
            return input_smiles, {}

        # Load the QM-predicted SDF
        if not os.path.exists(out_sdf):
            print(f"[{name}] ERROR: QupKake output file not found at {out_sdf}")
            print(f"[{name}] Falling back to input SMILES")
            return input_smiles, {}

        supplier = Chem.SDMolSupplier(out_sdf)
        mol = next(supplier, None)
        if mol is None:
            print(f"[{name}] ERROR: QupKake failed to generate a valid SDF")
            print(f"[{name}] Falling back to input SMILES")
            return input_smiles, {}

        # Extract pKa values from molecule-level properties (QupKake stores them there)
        pka_records = {}
        charge_adjustments = []

        mol_props = mol.GetPropsAsDict()

        # Get atom index and pKa value from molecule properties
        atom_idx_raw = mol_props.get('idx', None)
        pka_raw = mol_props.get('pka', None)
        pka_type = mol_props.get('pka_type', '')

        if atom_idx_raw is not None and pka_raw is not None:
            try:
                # Convert to strings for processing
                atom_idx_str = str(atom_idx_raw).strip()
                pka_str = str(pka_raw).strip()

                if not pka_str or pka_str.lower() in ['none', 'nan', '']:
                    pass
                else:
                    # Handle tensor format: "tensor(3.6062)" -> 3.6062
                    if 'tensor(' in pka_str:
                        pka_str = pka_str.split('(')[1].split(')')[0]

                    pka_val = float(pka_str)
                    atom_idx = int(atom_idx_str)

                # Get the actual atom
                if atom_idx < mol.GetNumAtoms():
                    atom = mol.GetAtomWithIdx(atom_idx)
                    atom_symbol = atom.GetSymbol()
                    original_charge = atom.GetFormalCharge()

                    # Store pKa value
                    pka_records[f"Atom {atom_idx} ({atom_symbol})"] = pka_val

                    # Basic Nitrogen (e.g., amines): Protonate if pKa > target pH
                    if atom.GetAtomicNum() == 7 and pka_val > target_ph:
                        atom.SetFormalCharge(1)
                        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                        if original_charge != 1:
                            charge_adjustments.append(f"N{atom_idx} (pKa={pka_val:.2f}): {original_charge:+d}→+1")

                    # Acidic Oxygen (e.g., carboxylic acids, ferulates): Deprotonate if pKa < target pH
                    elif atom.GetAtomicNum() == 8 and pka_val < target_ph:
                        atom.SetFormalCharge(-1)
                        atom.SetNumExplicitHs(max(0, atom.GetNumExplicitHs() - 1))
                        if original_charge != -1:
                            charge_adjustments.append(f"O{atom_idx} (pKa={pka_val:.2f}): {original_charge:+d}→-1")

            except (ValueError, IndexError):
                pass

        # Update valency and return the final physical SMILES
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            print(f"[{name}] ERROR: Sanitization failed after protonation adjustment: {e}")
            print(f"[{name}] Falling back to input SMILES")
            return input_smiles, {}

        corrected_smiles = Chem.MolToSmiles(mol)

        # Calculate and report charge changes
        mol_input = Chem.MolFromSmiles(input_smiles)
        mol_corrected = Chem.MolFromSmiles(corrected_smiles)

        if mol_input and mol_corrected:
            charge_input = Chem.GetFormalCharge(mol_input)
            charge_corrected = Chem.GetFormalCharge(mol_corrected)

            if charge_adjustments:
                print(f"[{name}] QupKake adjusted charge: {charge_input:+d} → {charge_corrected:+d}")
                for adj in charge_adjustments:
                    print(f"  • {adj}")
            else:
                print(f"[{name}] QupKake confirmed charge: {charge_corrected:+d}")

        print(f"[{name}] QupKake-corrected SMILES: {corrected_smiles}")
        return corrected_smiles, pka_records


# ============================================================================
# DEPRECATED: GNN-based pKa prediction (original implementation)
# This code has been moved here from qm_ligand_prep.py for reference
# ============================================================================

# try:
#     from pka_predictor.predict import predict
#     GNN_PKA_AVAILABLE = True
# except ImportError as e:
#     GNN_PKA_AVAILABLE = False
#     print(f"WARNING: pka_predictor not installed. Using input SMILES as-is. Error: {e}")
#
# def get_correct_protonation_state_gnn(input_smiles: str, target_ph: float = 7.4, name: str = "ligand") -> str:
#     """
#     Uses GNN pKa predictor to determine the dominant microstate at target pH.
#     Falls back to input SMILES if GNN is unavailable.
#
#     Args:
#         input_smiles: Input SMILES (any protonation state)
#         target_ph: Target pH (default 7.4)
#         name: Ligand name for logging
#
#     Returns:
#         SMILES of dominant microstate at target pH
#     """
#     if not GNN_PKA_AVAILABLE:
#         print(f"[{name}] GNN pKa predictor unavailable, using input SMILES")
#         return input_smiles
#
#     print(f"[{name}] Analyzing protonation state via GNN at pH {target_ph}...")
#
#     try:
#         import sys
#         old_argv = sys.argv
#         sys.argv = [sys.argv[0]]
#         try:
#             # predict() returns a DataFrame with 'protonated_smiles' column
#             df = predict(input_smiles, pH=target_ph)
#         finally:
#             sys.argv = old_argv
#
#         if df is None or df.empty or 'protonated_smiles' not in df.columns:
#             print(f"[{name}] GNN prediction failed or returned no results, using input SMILES")
#             return input_smiles
#
#         corrected_smiles = df.iloc[0]['protonated_smiles']
#
#         # Calculate formal charge difference
#         mol_input = Chem.MolFromSmiles(input_smiles)
#         mol_corrected = Chem.MolFromSmiles(corrected_smiles)
#
#         if mol_input and mol_corrected:
#             charge_input = Chem.GetFormalCharge(mol_input)
#             charge_corrected = Chem.GetFormalCharge(mol_corrected)
#
#             if charge_input != charge_corrected:
#                 print(f"[{name}] GNN adjusted charge: {charge_input:+d} → {charge_corrected:+d}")
#             else:
#                 print(f"[{name}] GNN confirmed charge: {charge_corrected:+d}")
#
#         print(f"[{name}] GNN-corrected SMILES: {corrected_smiles}")
#         return corrected_smiles
#
#     except Exception as e:
#         print(f"[{name}] ERROR during GNN prediction: {e}")
#         print(f"[{name}] Falling back to input SMILES")
#         return input_smiles
