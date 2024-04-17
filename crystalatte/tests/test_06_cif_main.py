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

    crystalatte.cif_main(fNameIn="crystalatte/data/cif/Ammonia.cif",
                         fNameOut="sc.xyz",
                         Na=1,
                         Nb=1,
                         Nc=1,
                         monomer_cutoff=0.0,
                         nmer_cutoff=0.0,
                         make_rect_box=True)

    ref = """16
        Crystal created from CIF file. Box size:       5.130500000000       5.130500000000       5.130500000000
        N                1.484253650000       4.049503650000       3.646246350000
        N                4.049503650000       3.646246350000       1.484253650000
        N                3.646246350000       1.484253650000       4.049503650000
        N                1.080996350000       1.080996350000       1.080996350000
        H                3.237858550000       3.935606550000       1.970625050000
        H                1.194893450000       4.535875050000       4.457891450000
        H                0.672608550000       3.760143450000       3.159874950000
        H                1.970625050000       3.237858550000       3.935606550000
        H                4.535875050000       4.457891450000       1.194893450000
        H                3.760143450000       3.159874950000       0.672608550000
        H                3.159874950000       0.672608550000       3.760143450000
        H                3.935606550000       1.970625050000       3.237858550000
        H                4.457891450000       1.194893450000       4.535875050000
        H                0.594624950000       1.892641450000       1.370356550000
        H                1.370356550000       0.594624950000       1.892641450000
        H                1.892641450000       1.370356550000       0.594624950000""".splitlines()


    with open("sc.xyz", "r") as xyz:
        for i, line in enumerate(xyz):
            assert line.strip() == ref[i].strip()

    # Clean-up generated test files.
    subprocess.call(["rm", "sc.xyz"])
