"""
Unit and regression test for the MolAligner package.
"""

# Import package, test suite, and other packages as needed
import MolAligner
import pytest
import sys

def test_MolAligner_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "MolAligner" in sys.modules
