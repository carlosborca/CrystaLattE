"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest
import subprocess

def test_less_than_dimers():
    """Checks that the program prints an error message when an invalid
    input with an 'nmers_up_to' keyword that is less than 2 is
    provided."""

    with pytest.raises(SystemExit):
        crystalatte.main(
            cif_input="crystalatte/data/cif/Carbon_Dioxide.cif", 
            nmers_up_to=1, 
            cle_run_type=["test"]
            )

    # Clean-up generated test files
    subprocess.call(["rm", "sc.xyz"])
