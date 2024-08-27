import io
import os
import pickle
import sys
import unittest
import textwrap
import tempfile

from scoring_matrices import ScoringMatrix

from .. import Alignment, Sequence, GappedSequence
from . import fasta, data

try:
    try:
        from importlib.resources import files as resource_files
    except ImportError:
        from importlib_resources import files as resource_files
except ImportError:
    resource_files = None


class TestAlignment(unittest.TestCase):

    def test_init(self):
        g1 = GappedSequence(b"test1", b"MNG-EGPNFY")
        g2 = GappedSequence(b"test2", b"MNGTEGP-FY")
        ali = Alignment([g1, g2])
        self.assertEqual(len(ali), 2)
        self.assertEqual(ali[0].id, g1.id)
        self.assertEqual(ali[1].id, g2.id)

        with self.assertRaises(TypeError):
            ali = Alignment(b"failure")
        with self.assertRaises(TypeError):
            ali = Alignment([b"failure"])

    def test_init_error(self):
        g1 = GappedSequence(b"test1", b"MNG-EGPNFY")
        g2 = GappedSequence(b"test2", b"MNGTEGP-FY")
        g3 = GappedSequence(b"test2", b"---TEGP-F")
        with self.assertRaises(ValueError):
            ali = Alignment([g1, g2, g3])

    def test_copy(self):
        g1 = GappedSequence(b"test1", b"MNG-EGPNFY")
        g2 = GappedSequence(b"test2", b"MNGTEGP-FY")
        ali = Alignment([g1, g2])
        copy = ali.copy()
        self.assertEqual(len(ali), len(copy))
        for i in range(len(ali)):
            self.assertEqual(ali[i].id, copy[i].id)
            self.assertEqual(ali[i].sequence, copy[i].sequence)

    def test_getitem(self):
        g1 = GappedSequence(b"test1", b"MNG-EGPNFY")
        g2 = GappedSequence(b"test2", b"MNGTEGP-FY")
        ali = Alignment([g1, g2])
        self.assertEqual(ali[0].sequence, b"MNG-EGPNFY")
        self.assertEqual(ali[1].sequence, b"MNGTEGP-FY")

    def test_pickle(self):
        g1 = GappedSequence(b"test1", b"MNG-EGPNFY")
        g2 = GappedSequence(b"test2", b"MNGTEGP-FY")
        ali = Alignment([g1, g2])
        pickled = pickle.loads(pickle.dumps(ali))
        self.assertEqual(len(ali), len(pickled))
        for i in range(len(ali)):
            self.assertEqual(ali[i].id, pickled[i].id)
            self.assertEqual(ali[i].sequence, pickled[i].sequence)

    def test_slice(self):
        sequences = [ GappedSequence(f"test{i}".encode(), b"M") for i in range(10) ]
        ali = Alignment(sequences)

        a1 = ali[:3]
        self.assertEqual(len(a1), 3)
        self.assertEqual([x.id for x in a1], [b"test0", b"test1", b"test2"])

        a2 = ali[-3:]
        self.assertEqual(len(a2), 3)
        self.assertEqual([x.id for x in a2], [b"test7", b"test8", b"test9"])

        a3 = ali[::2]
        self.assertEqual(len(a3), 5)
        self.assertEqual([x.id for x in a3], [b"test0", b"test2", b"test4", b"test6", b"test8"])