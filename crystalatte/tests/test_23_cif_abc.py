"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest

def test_noninteger_cif_a_input():
    """Checks that the program prints an error message when an invalid
    output with a non-integer 'cif_a' keyword is provided."""

    with pytest.raises(SystemExit):
        with pytest.raises(IndexError):
            crystalatte.input_parser("crystalatte/data/cle/cif_abc.cle")
