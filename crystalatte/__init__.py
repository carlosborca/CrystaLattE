"""
CrystaLattE
Automated calculation of crystal lattice energies with the many-body expansion
"""

# CrystaLattE imports
__all__ = ["input_parser",
        "extract_element",
        "write_xyz",
        "read_cif",
        "cif_main",
        "cif_driver",
        "center_supercell",
        "supercell2monomers",
        "create_nmer",
        "center_of_mass",
        "distance_matrix",
        "nre",
        "chemical_space",
        "build_nmer",
        "nmer2psiapimol",
        "monomer2makefp",
        "nmer2psithon",
        "psi4api_energies",
        "cle_manager",
        "print_results",
        "print_end_msg",
        "main"]

from .crystalatte import *

# Psithonyzer imports
__all__ = ["psithonyzer_success_check",
        "psithonyzer_get_nmer_data",
        "psithonyzer_print_header",
        "psithonyzer_print_results",
        "psithonyzer_print_end_msg",
        "psithonyzer_main"]

from .psithonyzer import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
