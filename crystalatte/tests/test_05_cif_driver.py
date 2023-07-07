"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest

def test_cle_cif_driver():
    """Checks that the CIF driver produces a correct list of arguments
    to be passed to the main CIF converter function."""

    ref = ['', '-i', 'crystalatte/data/cif/Ammonia.cif', '-o', 'sc.xyz', '-b', '3', '3', '3', '-r', '-d_mono', '0', '-d_nmer', '0']
    args = crystalatte.cif_driver("crystalatte/data/cif/Ammonia.cif", "sc.xyz", "3", "3", "3", "0", "0", 2)

    assert compare(ref, args)
