"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
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

    n_mono = len([k for k in nmers_dictionary.keys() if k.startswith("1mer-")])
    n_di = len([k for k in nmers_dictionary.keys() if k.startswith("2mer-")])
    n_tri = len([k for k in nmers_dictionary.keys() if k.startswith("3mer-")])
    n_tetra = len([k for k in nmers_dictionary.keys() if k.startswith("4mer-")])
    n_penta = len([k for k in nmers_dictionary.keys() if k.startswith("5mer-")])
    
    assert compare(7, n_mono, "Number of Monomers: ")
    assert compare(1, n_di, "Number of Monomers: ")
    assert compare(2, n_tri, "Number of Monomers: ")
    assert compare(1, n_tetra, "Number of Monomers: ")
    assert compare(1, n_penta, "Number of Monomers: ")

    for k, v in nmers_dictionary.items():
        n = int(k[0])
        ref = [4]*n
        if n != 1:
            assert compare(ref, v["atoms_per_monomer"])

    assert compare_values(-3.28520798120735e-07, nmers_dictionary["5mer-0+1+2+3+4"]["nambe"], atol=1.e-8)

    assert 0

