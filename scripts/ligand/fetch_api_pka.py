#!/usr/bin/env python3
"""
Fetch ChemAxon pKa data from ChEMBL API
Uses SMILES → InChIKey to query ChEMBL for experimental ChemAxon pKa values.
"""
import argparse
import sys
import requests
from rdkit import Chem

def get_chembl_pka(smiles, target_ph):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("ERROR: Invalid SMILES string.")
        sys.exit(1)

    inchikey = Chem.MolToInchiKey(mol)

    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_structures__standard_inchi_key={inchikey}"

    try:
        response = requests.get(url, headers={"Accept": "application/json"}, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f"ERROR: API request failed. {e}")
        sys.exit(1)

    data = response.json()
    if not data.get("molecules"):
        print(f"ERROR: Molecule not found in ChEMBL (InChIKey: {inchikey})")
        sys.exit(1)

    molecule_data = data["molecules"][0]
    chembl_id = molecule_data.get("molecule_chembl_id", "Unknown")
    pref_name = molecule_data.get("pref_name", "None")
    props = molecule_data.get("molecule_properties", {})

    if not props:
        print(f"ERROR: ChEMBL ID {chembl_id} found, but no property data available.")
        sys.exit(1)

    bpka_raw = props.get("cx_most_bpka")
    apka_raw = props.get("cx_most_apka")

    bpka = float(bpka_raw) if bpka_raw is not None else None
    apka = float(apka_raw) if apka_raw is not None else None

    charge = 0
    if bpka and bpka > target_ph:
        charge += 1
    if apka and apka < target_ph:
        charge -= 1

    print("-" * 40)
    print(f"SMILES:     {smiles}")
    print(f"InChIKey:   {inchikey}")
    print(f"Name:       {pref_name}")
    print(f"ChEMBL ID:  {chembl_id}")
    print("-" * 40)
    print(f"Target pH:            {target_ph}")
    print(f"Strongest Basic pKa:  {bpka if bpka is not None else 'None'}")
    print(f"Strongest Acidic pKa: {apka if apka is not None else 'None'}")
    print(f"Net Formal Charge:    {charge:+d}")
    print("-" * 40)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fetch exact ChemAxon pKa data from ChEMBL via SMILES.")
    parser.add_argument("smiles", help="Input SMILES string")
    parser.add_argument("--ph", type=float, default=7.4, help="Target pH (default: 7.4)")
    args = parser.parse_args()

    get_chembl_pka(args.smiles, args.ph)
