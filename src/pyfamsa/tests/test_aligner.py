import io
import os
import sys
import unittest
import textwrap
import tempfile

from scoring_matrices import ScoringMatrix

from .. import Aligner, Alignment, GappedSequence, Sequence
from . import fasta, data

try:
    try:
        from importlib.resources import files as resource_files
    except ImportError:
        from importlib_resources import files as resource_files
except ImportError:
    resource_files = None


class _Test(object):
     
    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_hemopexin_medoid_nj(self):
        self._test_famsa("hemopexin", "nj", "medoid")

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_hemopexin_medoid_sl(self):
        self._test_famsa("hemopexin", "sl", "medoid")

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_hemopexin_medoid_upgma(self):
        self._test_famsa("hemopexin", "upgma", "medoid")

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_adeno_fiber_medoid_upgma(self):
        self._test_famsa("adeno_fiber", "upgma", None)

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_adeno_fiber_sl(self):
        self._test_famsa("adeno_fiber", "sl", None)

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
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
        with resource_files(data).joinpath(filename).open() as file:
            records = list(fasta.parse(file))

        if tree_heuristic is None:
            filename = "{}.{}.afa".format(test_case, guide_tree)
        else:
            filename = "{}.{}-{}.afa".format(test_case, tree_heuristic, guide_tree)
        with resource_files(data).joinpath(filename).open() as file:
            result = list(fasta.parse(file))

        aligner = Aligner(guide_tree=guide_tree, tree_heuristic=tree_heuristic, threads=1)
        sequences = (
            Sequence(record.id.encode(), record.seq.encode())
            for record in records
        )

        alignment = aligner.align(sequences)
        for expected, actual in zip(result, alignment):
            self.assertEqual(expected.id, actual.id.decode())
            self.assertEqual(expected.seq, actual.sequence.decode())


class TestAlignProfiles(unittest.TestCase):

    @unittest.skipUnless(resource_files, "tests require `importlib.resources.files`")
    def test_adeno_fiber_upgma(self):
        with resource_files(data).joinpath("adeno_fiber.p1.afa").open() as file:
            a1 = Alignment(
                GappedSequence(record.id.encode(), record.seq.encode())
                for record in fasta.parse(file)
            )
        with resource_files(data).joinpath("adeno_fiber.p2.afa").open() as file:
            a2 = Alignment(
                GappedSequence(record.id.encode(), record.seq.encode())
                for record in fasta.parse(file)
            )

        aligner = Aligner(guide_tree="upgma", refine=False)
        alignment = aligner.align_profiles(a2, a1)

        with resource_files(data).joinpath("adeno_fiber.upgma.pp.afa").open() as file:
            result = list(fasta.parse(file))

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
        with resource_files(data).joinpath(filename).open() as file:
            records = list(fasta.parse(file))

        if tree_heuristic is None:
            filename = "{}.{}.nwk".format(test_case, guide_tree)
        else:
            filename = "{}.{}-{}.nwk".format(test_case, tree_heuristic, guide_tree)
        with resource_files(data).joinpath(filename).open() as file:
            result = file.read()

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
