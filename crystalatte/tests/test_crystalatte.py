"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest
import sys

def test_crystalatte_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "crystalatte" in sys.modules

def test_asdf():
    import crystalatte
    crystalatte.main(cif_input="../../Tests/Benzene/Benzene.cif", cif_a=2, cif_b=2, cif_c=2, nmers_up_to=2, r_cut_com=5.1, r_cut_monomer=3.3, r_cut_dimer=2.7, r_cut_trimer=2.7, r_cut_tetramer=2.7, r_cut_pentamer=5.6, cle_run_type=["psi4api", "quiet"], psi4_method="HF/STO-3G", psi4_bsse="cp", psi4_memory="500 MB", verbose=2)
