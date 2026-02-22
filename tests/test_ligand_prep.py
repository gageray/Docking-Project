import os
import pytest
from rdkit import Chem
import dimorphite_dl
import sys

# Add scripts directory to path to import ligand_prep
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'scripts'))
from ligand_prep import protonate, identify_tail_atoms, set_dihedral

def test_protonate_basic():
    """Test that DimorphiteDL protonation works on a simple amine at pH 7.4."""
    # Propylamine SMILES: CCCN
    smiles = "CCCN"
    result = protonate(smiles)
    
    # At pH 7.4, the primary amine should be protonated: [NH3+]CCC
    mol = Chem.MolFromSmiles(result)
    assert mol is not None, "Protonated SMILES should be valid"
    
    # Check formal charge is +1
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    assert charge == 1, f"Expected charge +1, got {charge}"

def test_protonate_acid():
    """Test protonation on a carboxylic acid at pH 7.4."""
    # Acetic acid SMILES: CC(=O)O
    smiles = "CC(=O)O"
    result = protonate(smiles)
    
    # At pH 7.4, it should be deprotonated: CC(=O)[O-]
    mol = Chem.MolFromSmiles(result)
    assert mol is not None
    
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    assert charge == -1, f"Expected charge -1, got {charge}"

def test_identify_tail_atoms():
    """Test the DFS tail atom identifier logic on a dummy ester."""
    # Dummy ester with a 4-carbon tail: CC(=O)OCCCC
    mol = Chem.MolFromSmiles("CC(=O)OCCCC")
    
    # SMARTS match for the ester linkage: [CX3H1]=[CX3H1]C(=O)[OX2][CH2]
    # Actually our script uses a specific match pattern in the code.
    # Let's manually define the ester match tuple for this test based on atomic indices.
    # C(0) - C(1)(=O(2)) - O(3) - C(4) - C(5) - C(6) - C(7)
    
    c_carbonyl = 1
    o_carbonyl = 2
    o_ester = 3
    c_alkyl1 = 4
    
    ester_match = (None, None, c_carbonyl, o_carbonyl, o_ester, c_alkyl1)
    
    tail_atoms = identify_tail_atoms(mol, ester_match)
    
    # The tail should include C(4), C(5), C(6), C(7).
    assert len(tail_atoms) == 4, f"Expected 4 tail carbons, found {len(tail_atoms)}"
    assert set(tail_atoms) == {4, 5, 6, 7}, f"Tail atom indices incorrect: {tail_atoms}"
