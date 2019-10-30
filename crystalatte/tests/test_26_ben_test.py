"""
Unit and regression test for the crystalatte package.
"""

# Import package, test suite, and other packages as needed
from qcelemental.testing import compare, compare_values
import crystalatte
import pytest
import sys

def test_testmode_benzene():
    """Test to reproduce the structures of A. L. Ringer and C. D. 
    Sherrill, Chem. Eur. J., 2008, 14, pp 2542–2547."""

    # Execute the main function of crystalatte and retrieve the N-mers dictionary.
    nmers, cle = crystalatte.main(cif_input="crystalatte/data/Benzene.cif", 
            cif_output="crystalatte/data/Benzene.xyz", 
            cif_a=5, 
            cif_b=5, 
            cif_c=5, 
            nmers_up_to=2, 
            r_cut_com=9.5, 
            r_cut_monomer=11.4, 
            r_cut_dimer=11.4, 
            r_cut_trimer=11.4, 
            r_cut_tetramer=11.4, 
            r_cut_pentamer=11.4, 
            cle_run_type=["test"], 
            psi4_method="HF/STO-3G", 
            psi4_bsse="nocp", 
            psi4_memory="500 MB", 
            verbose=2)
    
    # For debug.
    #import pprint
    #pprint.pprint(nmers)

    # Test the number of N-mers.
    number_mono  = len([k for k in nmers.keys() if k.startswith("1mer-")])
    number_di    = len([k for k in nmers.keys() if k.startswith("2mer-")])
    
    assert compare(84, number_mono,  "Number of Monomers: ")
    assert compare(10, number_di,    "Number of Monomers: ")

    # Test the number of atoms per monomer in each N-mer.
    for k, v in nmers.items():
        n = int(k[0])
        ref = [12]*n
        if n != 1:
            assert compare(ref, v["atoms_per_monomer"])

    # Test replicas for each N-mer.
    assert compare_values(4, nmers["2mer-0+1"]["replicas"],  "2mer-0+1 Replicas: ")
    assert compare_values(4, nmers["2mer-0+2"]["replicas"],  "2mer-0+1 Replicas: ")
    assert compare_values(4, nmers["2mer-0+6"]["replicas"],  "2mer-0+1 Replicas: ")
    assert compare_values(2, nmers["2mer-0+16"]["replicas"], "2mer-0+1 Replicas: ")
    assert compare_values(2, nmers["2mer-0+5"]["replicas"],  "2mer-0+1 Replicas: ")
    assert compare_values(2, nmers["2mer-0+21"]["replicas"], "2mer-0+1 Replicas: ")
    assert compare_values(4, nmers["2mer-0+12"]["replicas"], "2mer-0+1 Replicas: ")
    assert compare_values(4, nmers["2mer-0+14"]["replicas"], "2mer-0+1 Replicas: ")
    assert compare_values(4, nmers["2mer-0+15"]["replicas"], "2mer-0+1 Replicas: ")
    assert compare_values(4, nmers["2mer-0+27"]["replicas"], "2mer-0+1 Replicas: ")

    # Test the non-additive many-body energies.
    assert compare_values(0.0, nmers["2mer-0+1"]["nambe"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+2"]["nambe"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+6"]["nambe"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+16"]["nambe"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+5"]["nambe"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+21"]["nambe"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+12"]["nambe"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+14"]["nambe"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+15"]["nambe"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+27"]["nambe"], atol=1.e-8)

    # Test contributions for each N-mer.
    assert compare_values(0.0, nmers["2mer-0+1"]["contrib"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+2"]["contrib"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+6"]["contrib"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+16"]["contrib"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+5"]["contrib"],  atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+21"]["contrib"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+12"]["contrib"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+14"]["contrib"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+15"]["contrib"], atol=1.e-8)
    assert compare_values(0.0, nmers["2mer-0+27"]["contrib"], atol=1.e-8)

    # Test the nuclear repulsion energies for each N-mer.
    assert compare_values(564.7978818358, nmers["2mer-0+1"]["nre"],  atol=1.e-8)
    assert compare_values(569.6323170325, nmers["2mer-0+2"]["nre"],  atol=1.e-8)
    assert compare_values(588.8876674569, nmers["2mer-0+6"]["nre"],  atol=1.e-8)
    assert compare_values(542.3197525768, nmers["2mer-0+16"]["nre"], atol=1.e-8)
    assert compare_values(507.4201406916, nmers["2mer-0+5"]["nre"],  atol=1.e-8)
    assert compare_values(532.0439341175, nmers["2mer-0+21"]["nre"], atol=1.e-8)
    assert compare_values(507.4677080555, nmers["2mer-0+12"]["nre"], atol=1.e-8)
    assert compare_values(510.6225302261, nmers["2mer-0+14"]["nre"], atol=1.e-8)
    assert compare_values(509.3453310551, nmers["2mer-0+15"]["nre"], atol=1.e-8)
    assert compare_values(505.1441258078, nmers["2mer-0+27"]["nre"], atol=1.e-8)

    # Test the COM-based priority for each N-mer.
    assert compare_values(0.0006907226, nmers["2mer-0+1"]["priority_com"],  atol=1.e-8)
    assert compare_values(0.0007548316, nmers["2mer-0+2"]["priority_com"],  atol=1.e-8)
    assert compare_values(0.0011681199, nmers["2mer-0+6"]["priority_com"],  atol=1.e-8)
    assert compare_values(0.0004692044, nmers["2mer-0+16"]["priority_com"], atol=1.e-8)
    assert compare_values(0.0001772763, nmers["2mer-0+5"]["priority_com"],  atol=1.e-8)
    assert compare_values(0.0003671724, nmers["2mer-0+21"]["priority_com"], atol=1.e-8)
    assert compare_values(0.0001783188, nmers["2mer-0+12"]["priority_com"], atol=1.e-8)
    assert compare_values(0.0001987883, nmers["2mer-0+14"]["priority_com"], atol=1.e-8)
    assert compare_values(0.0001987883, nmers["2mer-0+15"]["priority_com"], atol=1.e-8)
    assert compare_values(0.0001783188, nmers["2mer-0+27"]["priority_com"], atol=1.e-8)

    # Test the N-mer-cutoff-based priority for each N-mer.
    assert compare_values(0.0090014493, nmers["2mer-0+1"]["priority_min"],  atol=1.e-8)
    assert compare_values(0.0081275168, nmers["2mer-0+2"]["priority_min"],  atol=1.e-8)
    assert compare_values(0.0077024945, nmers["2mer-0+6"]["priority_min"],  atol=1.e-8)
    assert compare_values(0.0014125011, nmers["2mer-0+16"]["priority_min"], atol=1.e-8)
    assert compare_values(0.0012169745, nmers["2mer-0+5"]["priority_min"],  atol=1.e-8)
    assert compare_values(0.0011442584, nmers["2mer-0+21"]["priority_min"], atol=1.e-8)
    assert compare_values(0.0010809527, nmers["2mer-0+12"]["priority_min"], atol=1.e-8)
    assert compare_values(0.0009506992, nmers["2mer-0+14"]["priority_min"], atol=1.e-8)
    assert compare_values(0.0006261051, nmers["2mer-0+15"]["priority_min"], atol=1.e-8)
    assert compare_values(0.0003528237, nmers["2mer-0+27"]["priority_min"], atol=1.e-8)

    # Test the N-mer cutoffs for each N-mer.
    assert compare_values(4.8072405390,  min(nmers["2mer-0+1"]["min_monomer_separations"]),  atol=1.e-8)
    assert compare_values(4.9737128407,  min(nmers["2mer-0+2"]["min_monomer_separations"]),  atol=1.e-8)
    assert compare_values(5.0635628061,  min(nmers["2mer-0+6"]["min_monomer_separations"]),  atol=1.e-8)
    assert compare_values(8.9125861155,  min(nmers["2mer-0+16"]["min_monomer_separations"]), atol=1.e-8)
    assert compare_values(9.3664029960,  min(nmers["2mer-0+5"]["min_monomer_separations"]),  atol=1.e-8)
    assert compare_values(9.5607501286,  min(nmers["2mer-0+21"]["min_monomer_separations"]), atol=1.e-8)
    assert compare_values(9.7438615012,  min(nmers["2mer-0+12"]["min_monomer_separations"]), atol=1.e-8)
    assert compare_values(10.1699533857, min(nmers["2mer-0+14"]["min_monomer_separations"]), atol=1.e-8)
    assert compare_values(11.6891856781, min(nmers["2mer-0+15"]["min_monomer_separations"]), atol=1.e-8)
    assert compare_values(14.1518785972, min(nmers["2mer-0+27"]["min_monomer_separations"]), atol=1.e-8)


    # Test the COM cutoffs for each N-mer.
    assert compare_values(11.3126785806, min(nmers["2mer-0+1"]["com_monomer_separations"]),  atol=1.e-8)
    assert compare_values(10.9828900652, min(nmers["2mer-0+2"]["com_monomer_separations"]),  atol=1.e-8)
    assert compare_values(9.4952015406,  min(nmers["2mer-0+6"]["com_monomer_separations"]),  atol=1.e-8)
    assert compare_values(12.8690349144, min(nmers["2mer-0+16"]["com_monomer_separations"]), atol=1.e-8)
    assert compare_values(17.8012201018, min(nmers["2mer-0+5"]["com_monomer_separations"]),  atol=1.e-8)
    assert compare_values(13.9650760671, min(nmers["2mer-0+21"]["com_monomer_separations"]), atol=1.e-8)
    assert compare_values(17.7664634564, min(nmers["2mer-0+12"]["com_monomer_separations"]), atol=1.e-8)
    assert compare_values(17.1344318930, min(nmers["2mer-0+14"]["com_monomer_separations"]), atol=1.e-8)
    assert compare_values(17.1344318930, min(nmers["2mer-0+15"]["com_monomer_separations"]), atol=1.e-8)
    assert compare_values(17.7664634564, min(nmers["2mer-0+27"]["com_monomer_separations"]), atol=1.e-8)

    # Test the crystal lattice energy.
    assert compare_values(0.0, cle, atol=1.e-8)
