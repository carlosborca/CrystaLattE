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

    l0 = None 
    l1 = None 
    l2 = None 
    l3 = None 
    l4 = None 
    l5 = None 
    l6 = None

    with open("sc.xyz", "r") as xyz:
        
        for l in xyz:
            
            if "5" in l:
                l0 = True

            if "Crystal created from CIF file. Box size:    1.00000    1.00000    1.00000" in l:
                l1 = True
            
            if "C            0.000000   0.000000   0.000000" in l:
                l2 = True

            if "H            0.640513  -0.640513   0.640513" in l:
                l3 = True

            if "H            0.640513   0.640513  -0.640513" in l:
                l4 = True

            if "H           -0.640513   0.640513   0.640513" in l:
                l5 = True

            if "H           -0.640513  -0.640513  -0.640513" in l:
                l6 = True

    assert compare(True, l0)
    assert compare(True, l1)
    assert compare(True, l2)
    assert compare(True, l3)
    assert compare(True, l4)
    assert compare(True, l5)
    assert compare(True, l6)

    # Clean-up generated test files.
    subprocess.call(["rm", "sc.xyz"])
