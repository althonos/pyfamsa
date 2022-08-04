# noqa: D104
from ._version import __version__  # isort: skip

from . import _famsa
from ._famsa import (
    Alignment,
    Aligner,
    Sequence,
    GappedSequence,
    Tree,
)

__doc__ = _famsa.__doc__
__all__ = [
    "Alignment",
    "Aligner",
    "Sequence",
    "GappedSequence"
]

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "GPLv3"
