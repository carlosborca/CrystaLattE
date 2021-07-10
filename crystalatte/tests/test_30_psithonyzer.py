"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import os
import pytest
import subprocess

def test_psithonyzer():
    """Checks that the CSV file produced by the psithonyzer script
    contains correct information."""

    # Get the root directory and chenge to the output directory.
    root = os.getcwd()
    d = os.path.join(root, "crystalatte", "data", "out")
    os.chdir(d)

    # Execute the psithonyzer script.
    subprocess.call(["./../../psithonyzer.py"])

    # Change directory back to root; do this before test,
    # or else a failed test may leave us in the wrong directory.
    os.chdir(root)

    csv_lines = []

    with open("crystalatte/data/out/Ammonia.csv", 'r') as csv:

        for line in csv:

            if "Energy" in line:
                continue
            
            else:
                csv_lines.append(line[:-1].replace(" ", ""))

    results = []

    results.append('2mer-0+1,0.87719178,6,2.63157535,2.63157535,8.552209e-03,3.8830,2.5880,2.588')
    results.append('3mer-0+1+2,0.03591565,3,0.03591565,2.66749100,7.570855e-07,3.7280,2.5347,2.428,2.588,2.588')
    results.append('3mer-0+1+5,-0.10631085,3,-0.10631085,2.56118015,6.255109e-07,3.8830,2.5880,2.588,2.588,2.588')
    results.append('4mer-0+1+2+3,0.00643113,1,0.00160778,2.56278793,6.937449e-13,3.6505,2.5080,2.428,2.428,2.428,2.588,2.588,2.588')
    results.append('5mer-0+1+2+3+4,-0.00086359,3,-0.00051816,2.56226977,1.441820e-22,4.1234,2.9181,2.428,2.428,2.428,2.588,2.588,2.588,2.588,2.588,3.945,5.012')

    assert compare(csv_lines, results)
    

    # Clean-up generated test files.
    subprocess.call(["rm", "crystalatte/data/out/Ammonia.csv"])
