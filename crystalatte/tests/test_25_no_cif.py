"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest

def test_no_cif_input():
    """Checks that the program prints an error message when an invalid
    input without a 'cif_input' keyword is provided."""

    with pytest.raises(SystemExit):
        with pytest.raises(IndexError):
            crystalatte.input_parser("crystalatte/data/cle/no_cif.cle")
