import datetime
import os
import typing

from scoring_matrices import ScoringMatrix

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

GuideTreeMethod = Literal["sl", "slink", "upgma", "nj"]
TreeHeuristicMethod = Literal["medoid", "part"]
Node = typing.Tuple[int, int]

__version__: str

FAMSA_ALPHABET: str
MIQS: ScoringMatrix

class _VersionInfo(typing.NamedTuple):
    major: int
    minor: int
    micro: int

class _Info(typing.NamedTuple):
    version: str
    version_info: _VersionInfo
    date: datetime.date
    authors: typing.List[str]

def famsa_info() -> _Info: ...

class Sequence:
    def __init__(self, id: bytes, sequence: bytes) -> None: ...
    def __copy__(self) -> Sequence: ...
    def __repr__(self) -> str: ...
    @property
    def id(self) -> bytes: ...
    @property
    def sequence(self) -> bytes: ...
    @property
    def size(self) -> int: ...
    def copy(self) -> Sequence: ...

class GappedSequence:
    def __init__(self, id: bytes, sequence: bytes) -> None: ...
    def __copy__(self) -> GappedSequence: ...
    def __repr__(self) -> str: ...
    @property
    def id(self) -> bytes: ...
    @property
    def sequence(self) -> bytes: ...
    @property
    def size(self) -> int: ...
    @property
    def gapped_size(self) -> int: ...
    def copy(self) -> GappedSequence: ...

class Alignment(typing.Sequence[GappedSequence]):
    def __init__(self, objects: typing.Iterable[GappedSequence] = ()) -> None: ...
    def __len__(self) -> int: ...
    @typing.overload
    def __getitem__(self, index: int) -> GappedSequence: ... 
    @typing.overload
    def __getitem__(self, index: slice) -> Alignment: ...
    @typing.overload
    def __getitem__(self, index: typing.Union[slice, int]) -> typing.Union[Alignment, GappedSequence]: ...

class Aligner:
    def __init__(
        self,
        *,
        threads: int = 0,
        guide_tree: GuideTreeMethod = "sl",
        tree_heuristic: typing.Optional[TreeHeuristicMethod] = None,
        medoid_threshold: int = 0,
        n_refinements: int = 100,
        keep_duplicates: bool = False,
        refine: typing.Optional[bool] = None,
        scoring_matrix: typing.Union[ScoringMatrix, str, None] = None,
    ) -> None: ...
    def align(self, sequences: typing.Iterable[Sequence]) -> Alignment: ...
    def align_profiles(self, profile1: Alignment, profile2: Alignment) -> Alignment: ...
    def build_tree(self, sequences: typing.Iterable[Sequence]) -> GuideTree: ...

class GuideTree(typing.Sequence[Node]):
    def __init__(self) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Node: ...  # type: ignore
    def dumps(self) -> bytes: ...
    def dump(
        self,
        file: typing.Union[
            str, bytes, os.PathLike[bytes], os.PathLike[str], typing.BinaryIO
        ],
    ) -> int: ...
