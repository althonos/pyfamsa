import os
import typing

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal  # type: ignore

GuideTreeMethod = Literal["sl", "slink", "upgma", "nj"]
TreeHeuristicMethod = Literal["medoid", "part"]
Node = typing.Tuple[int, int]

class Sequence:
    def __init__(self, id: bytes, sequence: bytes) -> None: ...
    def __copy__(self) -> Sequence: ...
    def __repr__(self) -> str: ...
    @property
    def id(self) -> bytes: ...
    @property
    def sequence(self) -> bytes: ...
    def copy(self) -> Sequence: ...

class GappedSequence:
    @property
    def id(self) -> bytes: ...
    @property
    def sequence(self) -> bytes: ...

class Alignment(typing.Sequence[GappedSequence]):
    def __init__(self) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> GappedSequence: ...  # type: ignore

class Aligner:
    def __init__(
        self,
        *,
        threads: int = 0,
        guide_tree: GuideTreeMethod = "sl",
        tree_heuristic: typing.Optional[TreeHeuristicMethod] = None,
        medoid_threshold: int = 0,
        n_refinements: int = 100,
        force_refinement: bool = False,
    ) -> None: ...
    def align(self, sequences: typing.Iterable[Sequence]) -> Alignment: ...
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
