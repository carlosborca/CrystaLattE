"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_psz_success_check():
    """Checks that the psithonyzer success checker function retruns a
    True boolean value when it finds a Psi4's successful execution
    string."""

    with open("beer.out", "w") as f:
        f.write("# This is a dummy Psi4 output test file.")
        f.write("Psi4 exiting successfully. Buy a developer a beer!")
    
    with open("coffee.out", "w") as f:
        f.write("Psi4 encountered an error. Buy a developer more coffee!")

    beer = crystalatte.psz_success_check("beer.out", 2)
    coffee = crystalatte.psz_success_check("coffee.out", 2)
    assert compare(True, beer)
    assert compare(False, coffee)

    # Clean-up generated test files.
    subprocess.call(["rm", "beer.out"])
    subprocess.call(["rm", "coffee.out"])
