"""
FAST: Scalable Similarity-based Test Case Prioritization

A Python library for test case prioritization using similarity-based algorithms
including LSH-based approaches and competitor algorithms.
"""

__version__ = "1.3.0"
__author__ = "Breno Miranda, Emilio Cruciani, Roberto Verdecchia, Antonia Bertolino"
__email__ = ""
__description__ = (
    "FAST Approaches to Scalable Similarity-based Test Case Prioritization"
)

from .prioritize import bbox_prioritization, run_blackbox_file

__all__ = ["bbox_prioritization", "run_blackbox_file"]
