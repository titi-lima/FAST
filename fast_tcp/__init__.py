"""
FAST: Scalable Similarity-based Test Case Prioritization

A Python library for test case prioritization using similarity-based algorithms
including LSH-based approaches and competitor algorithms.
"""

__version__ = "1.0.0"
__author__ = "Breno Miranda, Emilio Cruciani, Roberto Verdecchia, Antonia Bertolino"
__email__ = ""
__description__ = (
    "FAST Approaches to Scalable Similarity-based Test Case Prioritization"
)

# Import main functionality
from .fast import fast_pw, fast_
from .competitors import str_, i_tsd, gt, ga, ga_s, artd, artf
from .metric import apfd
from .lsh import kShingles, tcMinhashing, LSHBucket, LSHCandidates, jDistanceEstimate

__all__ = [
    "fast_pw",
    "fast_",
    "str_",
    "i_tsd",
    "gt",
    "ga",
    "ga_s",
    "artd",
    "artf",
    "apfd",
    "kShingles",
    "tcMinhashing",
    "LSHBucket",
    "LSHCandidates",
    "jDistanceEstimate",
]
