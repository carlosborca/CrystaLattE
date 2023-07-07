"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

@pytest.mark.skip(reason="This error was downgraded to a warning.")
def test_monomer_cutoff_co2():
    """Checks that the program prints an error message when an invalid
    monomer cutoff, that is loger than half the smallest dimension of
    the supercell has been chosen."""

    with pytest.raises(SystemExit):
        crystalatte.main(
            cif_input="crystalatte/data/cif/Carbon_Dioxide.cif", 
            cif_a=3, 
            cif_b=3, 
            cif_c=3, 
            nmers_up_to=2, 
            r_cut_monomer=8.2, 
            cle_run_type=["test"]
            )
    
    # Clean-up generated test files
    subprocess.call(["rm", "sc.xyz"])
