import io
import os
import sys
import unittest
import textwrap
import tempfile

from .. import Aligner, Sequence
from . import fasta, data

try:
    try:
        import importlib.resources as importlib_resources
    except ImportError:
        import importlib_resources
except ImportError:
    importlib_resources = None


class _Test(object):

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_hemopexin_medoid_nj(self):
        self._test_famsa("hemopexin", "nj")

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_hemopexin_medoid_sl(self):
        self._test_famsa("hemopexin", "sl")

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_hemopexin_medoid_upgma(self):
        self._test_famsa("hemopexin", "upgma")

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_adeno_fiber_medoid_sl(self):
        self._test_famsa("adeno_fiber", "sl")

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_adeno_fiber_medoid_upgma(self):
        self._test_famsa("adeno_fiber", "upgma")

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_adeno_fiber_sl(self):
        self._test_famsa("adeno_fiber", "sl", None)

    @unittest.skipUnless(importlib_resources, "tests require `importlib.resources`")
    def test_adeno_fiber_upgma(self):
        self._test_famsa("adeno_fiber", "upgma", None)


class TestAlign(unittest.TestCase, _Test):

    def _test_famsa(self, test_case, guide_tree, tree_heuristic="medoid"):
        filename = "{}.faa".format(test_case)
        with importlib_resources.open_text(data.__name__, filename) as file:
            records = list(fasta.parse(file))

        if tree_heuristic is None:
            filename = "{}.{}.afa".format(test_case, guide_tree)
        else:
            filename = "{}.{}-{}.afa".format(test_case, tree_heuristic, guide_tree)
        with importlib_resources.open_text(data.__name__, filename) as file:
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


class TestBuildTree(unittest.TestCase, _Test):

    def _test_famsa(self, test_case, guide_tree, tree_heuristic="medoid"):
        filename = "{}.faa".format(test_case)
        with importlib_resources.open_text(data.__name__, filename) as file:
            records = list(fasta.parse(file))

        if tree_heuristic is None:
            filename = "{}.{}.nwk".format(test_case, guide_tree)
        else:
            filename = "{}.{}-{}.nwk".format(test_case, tree_heuristic, guide_tree)
        with importlib_resources.open_text(data.__name__, filename) as file:
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
