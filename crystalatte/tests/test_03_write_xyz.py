"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_cle_write_xyz():
    """"Checks that the function to write xyz files from given atomic 
    coordinates, box size, and an output file name is able to correctly
    generate a new xyz file."""

    atoms = [("C", 0.0000000, 0.0000000, 0.0000000),
            ("H", 0.6405128, -0.6405128, 0.6405128),
            ("H", 0.6405128, 0.6405128, -0.6405128),
            ("H", -0.6405128, 0.6405128, 0.6405128),
            ("H", -0.6405128, -0.6405128, -0.6405128)]
    
    box = (1.0, 1.0, 1.0)

    with open("sc.xyz", "w") as f:
        crystalatte.write_xyz(atoms, box, f)

    ref = """5
        Crystal created from CIF file. Box size:       1.000000000000       1.000000000000       1.000000000000
        C                0.000000000000       0.000000000000       0.000000000000
        H                0.640512800000      -0.640512800000       0.640512800000
        H                0.640512800000       0.640512800000      -0.640512800000
        H               -0.640512800000       0.640512800000       0.640512800000
        H               -0.640512800000      -0.640512800000      -0.640512800000
        """.splitlines()

    with open("sc.xyz", "r") as xyz:
        for i, line in enumerate(xyz):
            assert line.strip() == ref[i].strip()

    # Clean-up generated test files.
    subprocess.call(["rm", "sc.xyz"])
