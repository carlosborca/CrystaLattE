"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_cle_cif_main():
    """Checks that the main CIF function can take give arguments and 
    produce a correct XYZ file containing a supercell based on the 
    provided CIF file for an ammonia crystal."""

    args = ['', '-i', 'crystalatte/data/cif/Ammonia.cif', '-o', 'sc.xyz', '-b', '1', '1', '1', '-r']

    crystalatte.cif_main(args)

    l00 = None
    l01 = None
    l02 = None
    l03 = None
    l04 = None
    l05 = None
    l06 = None
    l07 = None
    l08 = None
    l09 = None
    l10 = None
    l11 = None
    l12 = None
    l13 = None
    l14 = None
    l15 = None
    l16 = None

    with open("sc.xyz", "r") as f:
        for l in f:
            if "16" in l:
                l00 = True
            if "Crystal created from CIF file. Box size:    5.13050    5.13050    5.13050" in l:
                l01 = True
            if "N            1.484254   4.049504   3.646246" in l:
                l02 = True
            if "N            4.049504   3.646246   1.484254" in l:
                l03 = True
            if "N            3.646246   1.484254   4.049504" in l:
                l04 = True
            if "N            1.080996   1.080996   1.080996" in l:
                l05 = True
            if "H            3.237859   3.935607   1.970625" in l:
                l06 = True
            if "H            1.194893   4.535875   4.457891" in l:
                l07 = True
            if "H            1.970625   3.237859   3.935607" in l:
                l08 = True
            if "H            4.535875   4.457891   1.194893" in l:
                l09 = True
            if "H            3.760143   3.159875   0.672609" in l:
                l10 = True
            if "H            3.159875   0.672609   3.760143" in l:
                l11 = True
            if "H            3.935607   1.970625   3.237859" in l:
                l12 = True
            if "H            4.457891   1.194893   4.535875" in l:
                l13 = True
            if "H            0.594625   1.892641   1.370357" in l:
                l14 = True
            if "H            1.370357   0.594625   1.892641" in l:
                l15 = True
            if "H            1.892641   1.370357   0.594625" in l:
                l16 = True

    assert compare(True, l00)
    assert compare(True, l01)
    assert compare(True, l02)
    assert compare(True, l03)
    assert compare(True, l04)
    assert compare(True, l05)
    assert compare(True, l06)
    assert compare(True, l07)
    assert compare(True, l08)
    assert compare(True, l09)
    assert compare(True, l10)
    assert compare(True, l11)
    assert compare(True, l12)
    assert compare(True, l13)
    assert compare(True, l14)
    assert compare(True, l15)
    assert compare(True, l16)

    # Clean-up generated test files.
    subprocess.call(["rm", "sc.xyz"])
