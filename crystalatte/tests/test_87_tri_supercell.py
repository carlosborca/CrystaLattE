"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import subprocess

def test_supercell_triazine():
    """Checks that the program prints a warning message indicating that
    the supercell contains duplicate coordinates and that they will be
    skipped. Then, it should state the new number of unique coordinates
    in the supercell."""

    # Execute the main function of crystalatte and retrieve the N-mers dictionary.
    nmers, cle = crystalatte.main(
            cif_input="crystalatte/data/cif/Triazine.cif", 
            cif_a=3, 
            cif_b=3, 
            cif_c=3, 
            nmers_up_to=2, 
            r_cut_monomer=8.0, 
            r_cut_dimer=2.9, 
            cle_run_type=["test"], 
            verbose=2
            )
    
    # For debug.
    #import pprint
    #pprint.pprint(nmers)

    # Test the number of N-mers.
    number_mono  = len([k for k in nmers.keys() if k.startswith("1mer-")])
    number_di    = len([k for k in nmers.keys() if k.startswith("2mer-")])
    
    assert compare(39, number_mono,  "Number of Monomers: ")
    assert compare(1, number_di,    "Number of Dimers: ")

    # Test the number of atoms per monomer in each N-mer.
    for k, v in nmers.items():
        n = int(k[0])
        ref = [9]*n
        if n != 1:
            assert compare(ref, v["atoms_per_monomer"])

    # Test replicas for each N-mer.
    assert compare_values(6, nmers["2mer-0+1"]["replicas"],  "2mer-0+1 Replicas: ")
    #assert compare_values(1, nmers["2mer-0+8"]["replicas"],  "2mer-0+8 Replicas: ")
    #assert compare_values(1, nmers["2mer-0+10"]["replicas"], "2mer-0+10 Replicas: ")

    # Test the nuclear repulsion energies for each N-mer.
    assert compare_values(603.1477732010773, nmers["2mer-0+1"]["nre"],  atol=1.e-8)
    #assert compare_values(603.1478268883, nmers["2mer-0+8"]["nre"],  atol=1.e-8)
    #assert compare_values(603.1478008660, nmers["2mer-0+10"]["nre"], atol=1.e-8)

    # Test the COM-based priority for each N-mer.
    assert compare_values(0.00080001155, nmers["2mer-0+1"]["priority_com"],  atol=1.e-8)
    #assert compare_values(0.0008000115, nmers["2mer-0+8"]["priority_com"],  atol=1.e-8)
    #assert compare_values(0.0008000115, nmers["2mer-0+10"]["priority_com"], atol=1.e-8)

    # Test the N-mer-separation-based priority for each N-mer.
    assert compare_values(0.006754970684989542, nmers["2mer-0+1"]["priority_min"],  atol=1.e-8)
    #assert compare_values(0.0067549681, nmers["2mer-0+8"]["priority_min"],  atol=1.e-8)
    #assert compare_values(0.0067549700, nmers["2mer-0+10"]["priority_min"], atol=1.e-8)

    # Test the N-mer-cutoff-based priority for each N-mer.
    assert compare_values(0.0012769227508507793, nmers["2mer-0+1"]["priority_cutoff"],  atol=1.e-8)
    #assert compare_values(0.0012769220, nmers["2mer-0+8"]["priority_cutoff"],  atol=1.e-8)
    #assert compare_values(0.0012769225, nmers["2mer-0+10"]["priority_cutoff"], atol=1.e-8)

    # Test the N-mer cutoffs for each N-mer.
    assert compare_values(5.2900386342782975,  min(nmers["2mer-0+1"]["min_monomer_separations"]),  atol=1.e-8)
    #assert compare_values(5.2900393127,  min(nmers["2mer-0+8"]["min_monomer_separations"]),  atol=1.e-8)
    #assert compare_values(5.2900388235,  min(nmers["2mer-0+10"]["min_monomer_separations"]), atol=1.e-8)

    ## Test the COM cutoffs for each N-mer.
    assert compare_values(10.772121612292835, min(nmers["2mer-0+1"]["com_monomer_separations"]),  atol=1.e-8)
    #assert compare_values(10.7721212935, min(nmers["2mer-0+8"]["com_monomer_separations"]),  atol=1.e-8)
    #assert compare_values(10.7721220245, min(nmers["2mer-0+10"]["com_monomer_separations"]), atol=1.e-8)

    # Test the crystal lattice energy.
    assert compare_values(0.0, cle, atol=1.e-8)
    
    # Clean-up generated test files
    subprocess.call(["rm", "sc.xyz"])
