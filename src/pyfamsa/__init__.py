# noqa: D104
from . import _famsa
from ._famsa import (
    Aligner,
    Alignment,
    GappedSequence,
    GuideTree,
    Sequence,
    famsa_info,
    FAMSA_ALPHABET,
    MIQS,
    PFASUM31,
    PFASUM43,
    PFASUM60,
)

__doc__ = _famsa.__doc__
__version__ = _famsa.__version__
__all__ = [
    "Aligner",
    "Alignment",
    "GappedSequence",
    "GuideTree",
    "Sequence",
    "famsa_info",
    "FAMSA_ALPHABET",
    "MIQS"
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
