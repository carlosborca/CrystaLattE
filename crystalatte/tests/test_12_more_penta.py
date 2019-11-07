"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_more_than_pentamers():
    """Checks that the program prints an error message when an invalid
    input with an 'nmers_up_to' keyword that is more than 5 is
    provided."""

    with pytest.raises(SystemExit):
        crystalatte.main(
            cif_input="crystalatte/data/cif/Carbon_Dioxide.cif", 
            nmers_up_to=6, 
            r_cut_monomer=3.0, 
            cle_run_type=["test"]
            )

    # Clean-up generated test files
    subprocess.call(["rm", "sc.xyz"])
