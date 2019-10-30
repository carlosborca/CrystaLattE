"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import sys

def test_more_than_pentamers():
    """Main test of the psi4api mode with the ammonia crystal."""

    with pytest.raises(SystemExit):
        crystalatte.main(cif_input="crystalatte/data/Carbon_Dioxide.cif", 
            cif_output="crystalatte/data/Carbon_Dioxide.xyz", 
            cif_a=3, 
            cif_b=3, 
            cif_c=3, 
            nmers_up_to=6, 
            r_cut_com=1.0,
            r_cut_monomer=3.0, 
            cle_run_type=["test"]
            ) 