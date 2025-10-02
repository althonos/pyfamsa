import io
import os
import sys
import unittest
import textwrap
import tempfile

from scoring_matrices import ScoringMatrix

from .. import Aligner, Alignment, GappedSequence, Sequence
from . import fasta

try:
    try:
        from importlib.resources import files as resource_files
    except ImportError:
        from importlib_resources import files as resource_files
except ImportError:
    resource_files = None

def _load(filename, parse_fasta=True):
    if resource_files is None:
        raise unittest.SkipTest("tests require `importlib.resources.files`")
    try:
        path = resource_files("pyfamsa.tests.data").joinpath(filename)
    except ImportError:
        raise unittest.SkipTest(f"missing data module `pyfamsa.tests.data`")
    if not path.exists():
        raise unittest.SkipTest(f"missing data file {filename!r}")
    with path.open() as file:
        if parse_fasta:
            yield from fasta.parse(file)
        else:
            yield from file



class _Test(object):
     
    def test_hemopexin_medoid_nj(self):
        self._test_famsa("hemopexin", "nj", "medoid")

    def test_hemopexin_medoid_sl(self):
        self._test_famsa("hemopexin", "sl", "medoid")

    def test_hemopexin_medoid_upgma(self):
        self._test_famsa("hemopexin", "upgma", "medoid")

    def test_adeno_fiber_medoid_upgma(self):
        self._test_famsa("adeno_fiber", "upgma", None)

    def test_adeno_fiber_sl(self):
        self._test_famsa("adeno_fiber", "sl", None)

    def test_adeno_fiber_upgma(self):
        self._test_famsa("adeno_fiber", "upgma", None)


class TestAligner(unittest.TestCase):

    def test_init_scoring_matrix_str(self):
        matrix = ScoringMatrix.from_name("BLOSUM62")
        aligner = Aligner(scoring_matrix="BLOSUM62")
        self.assertEqual(aligner.scoring_matrix, matrix)

    def test_init_scoring_matrix_object(self):
        matrix = ScoringMatrix.from_name("BLOSUM62")
        aligner = Aligner(scoring_matrix=matrix)
        self.assertEqual(aligner.scoring_matrix, matrix)

    def test_init_scoring_matrix_error(self):
        self.assertRaises(TypeError, Aligner, scoring_matrix=1)


class TestAlign(unittest.TestCase, _Test):

    def test_no_sequence(self):
        aligner = Aligner()
        alignment = aligner.align([])
        self.assertEqual(len(alignment), 0)

    def test_single_sequence(self):
        seq = Sequence(b"test", b"ATGC")
        aligner = Aligner()
        alignment = aligner.align([seq])
        self.assertEqual(len(alignment), 1)
        self.assertEqual(alignment[0].sequence, b"ATGC")

    def _test_famsa(self, test_case, guide_tree, tree_heuristic):
        filename = "{}.faa".format(test_case)
        records = list(_load(filename))

        if tree_heuristic is None:
            filename = "{}.{}.afa".format(test_case, guide_tree)
        else:
            filename = "{}.{}-{}.afa".format(test_case, tree_heuristic, guide_tree)
        result = list(_load(filename))

        aligner = Aligner(guide_tree=guide_tree, tree_heuristic=tree_heuristic, threads=1)
        sequences = (
            Sequence(record.id.encode(), record.seq.encode())
            for record in records
        )

        alignment = aligner.align(sequences)
        for expected, actual in zip(result, alignment):
            self.assertEqual(expected.id, actual.id.decode())
            self.assertEqual(expected.seq, actual.sequence.decode())

    def test_scoring_matrix_smaller(self):

        scoring_matrix = ScoringMatrix.from_name("NUC.4.4").shuffle("ATGCN")
        aligner = Aligner(scoring_matrix=scoring_matrix)

        s1 = Sequence(b't1', b'ATGC')
        s2 = Sequence(b't2', b'ATTC')
        aligner.align([s1, s2])

    def test_scoring_matrix_smaller_invalid_character(self):

        scoring_matrix = ScoringMatrix.from_name("NUC.4.4").shuffle("ATGCN")
        aligner = Aligner(scoring_matrix=scoring_matrix)

        s1 = Sequence(b't1', b'ATGC')
        s2 = Sequence(b't2', b'ATTY')
        
        with self.assertRaises(ValueError) as ctx:
            aligner.align([s1, s2])


class TestAlignProfiles(unittest.TestCase):

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_adeno_fiber_upgma(self):
        a1 = Alignment(
            GappedSequence(record.id.encode(), record.seq.encode())
            for record in _load("adeno_fiber.p1.afa")
        )
        a2 = Alignment(
            GappedSequence(record.id.encode(), record.seq.encode())
            for record in _load("adeno_fiber.p2.afa")
        )

        aligner = Aligner(guide_tree="upgma", refine=False)
        alignment = aligner.align_profiles(a2, a1)

        result = list(_load("adeno_fiber.upgma.pp.afa"))
        for expected, actual in zip(result, alignment):
            self.assertEqual(expected.id, actual.id.decode())
            self.assertEqual(expected.seq, actual.sequence.decode())


class TestBuildTree(unittest.TestCase, _Test):

    def test_no_sequence(self):
        aligner = Aligner()
        tree = aligner.build_tree([])
        self.assertEqual(len(tree), 0)

    def test_single_sequence(self):
        seq = Sequence(b"test", b"ATGC")
        aligner = Aligner()
        tree = aligner.build_tree([seq])
        self.assertEqual(len(tree), 1)

    def _test_famsa(self, test_case, guide_tree, tree_heuristic):
        filename = "{}.faa".format(test_case)
        records = list(_load(filename))

        if tree_heuristic is None:
            filename = "{}.{}.nwk".format(test_case, guide_tree)
        else:
            filename = "{}.{}-{}.nwk".format(test_case, tree_heuristic, guide_tree)
        result = ''.join(_load(filename, parse_fasta=False))

        aligner = Aligner(guide_tree=guide_tree, tree_heuristic=tree_heuristic, threads=1)
        sequences = (
            Sequence(record.id.encode(), record.seq.encode())
            for record in records
        )

        tree = aligner.build_tree(sequences)
        self.assertMultiLineEqual(
            tree.dumps().decode(),
            result,
        )
