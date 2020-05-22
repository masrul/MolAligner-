"""
MolAligner 
Align a molecule along a desired direction
"""

# Add imports here
from .molaligner import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
