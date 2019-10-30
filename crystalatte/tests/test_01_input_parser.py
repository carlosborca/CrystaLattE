"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import sys

def test_input_parser():
    """
    ."""

    # Execute the main function of crystalatte and retrieve the N-mers dictionary.
    keywords = crystalatte.input_parser("crystalatte/data/Ammonia.cle")
    
    # For debug.
    #import pprint
    #pprint.pprint(keywords)

    assert compare_values(3, keywords['cif_a'])
    assert compare_values(3, keywords['cif_b'])
    assert compare_values(3, keywords['cif_c'])
    assert compare_values(1.2, keywords['bfs_thresh'])
    assert compare_values(2, keywords['nmers_up_to'])
    assert compare_values(6.5, keywords['r_cut_com'])
    assert compare_values(2.6, keywords['r_cut_dimer'])
    assert compare_values(3.5, keywords['r_cut_monomer'])
    assert compare_values(6.1, keywords['r_cut_pentamer'])
    assert compare_values(3.7, keywords['r_cut_tetramer'])
    assert compare_values(3.7, keywords['r_cut_trimer'])
    assert compare_values(2, keywords['verbose'])
    assert compare("crystalatte/data/Ammonia.cif", keywords['cif_input'])
    assert compare("crystalatte/data/Ammonia.xyz", keywords['cif_output'])
    assert compare(["test"], keywords['cle_run_type'])
    assert compare("ChSEV", keywords['uniq_filter'])
    assert compare("nocp", keywords['psi4_bsse'])
    assert compare("500 MB", keywords['psi4_memory'])
    assert compare('HF/STO-3G', keywords['psi4_method'])
