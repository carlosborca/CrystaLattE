"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import sys

def test_makefp_benzene():
    """
    ."""

    # At the moment this test only checks that makefp mode ends.
    with pytest.raises(SystemExit):
        crystalatte.main(cif_input="crystalatte/data/Benzene.cif", 
            cif_a=3, 
            cif_b=3, 
            cif_c=3, 
            nmers_up_to=2, 
            r_cut_monomer=2.0, 
            cle_run_type=["makefp"], 
            verbose=2)
