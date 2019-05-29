"""
CrystaLattE
Automated calculation of crystal lattice energies with the many-body expansion
"""

# Add imports here
from .crystalatte import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
