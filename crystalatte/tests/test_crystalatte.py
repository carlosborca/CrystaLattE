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
    import pprint
    nmers_dictionary = crystalatte.main(cif_input="../Tests/Ammonia/Ammonia.cif", cif_output="../Tests/Ammonia/Ammonia.xyz", cif_a=3, cif_b=3, cif_c=3, nmers_up_to=5, r_cut_com=6.5, r_cut_monomer=3.5, r_cut_dimer=2.6, r_cut_trimer=3.7, r_cut_tetramer=3.7, r_cut_pentamer=6.1, cle_run_type=["psi4api"], psi4_method="HF/STO-3G", psi4_bsse="nocp", psi4_memory="500 MB", verbose=2)
    pprint.pprint(nmers_dictionary)
