"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
import crystalatte
import pytest
import sys

def test_crystalatte_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "crystalatte" in sys.modules
