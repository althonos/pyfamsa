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


class TestAligner(unittest.TestCase):

    def _test_hemopexin(self, guide_tree, tree_heuristic="medoid"):
        with importlib_resources.open_text(data.__name__, "hemopexin.faa") as file:
            hemopexin = list(fasta.parse(file))
        with importlib_resources.open_text(data.__name__, "hemopexin.{}-{}.afa".format(tree_heuristic, guide_tree)) as file:
            result = list(fasta.parse(file))

        aligner = Aligner(guide_tree=guide_tree, tree_heuristic=tree_heuristic)
        sequences = (
            Sequence(record.id.encode(), record.seq.encode())
            for record in hemopexin
        )

        alignment = aligner.align(sequences)
        for expected, actual in zip(result, alignment):
            self.assertEqual(expected.id, actual.id.decode())
            self.assertEqual(expected.seq, actual.sequence.decode())

    def test_hemopexin_nj(self):
        self._test_hemopexin("nj")

    def test_hemopexin_sl(self):
        self._test_hemopexin("sl")

    def test_hemopexin_upgma(self):
        self._test_hemopexin("upgma")
