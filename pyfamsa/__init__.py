# noqa: D104
from ._version import __version__  # isort: skip

from . import _famsa
from ._famsa import (
    Aligner,
    Alignment,
    GappedSequence,
    GuideTree,
    Sequence,
    famsa_info,
)

__doc__ = _famsa.__doc__
__all__ = [
    "Aligner",
    "Alignment",
    "GappedSequence"
    "GuideTree"
    "Sequence",
    "famsa_info",
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
