import unittest

from .. import Sequence, GappedSequence


class TestSequence(unittest.TestCase):

    def test_sequence_property(self):
        seq = Sequence(b"test", b"MYYK")
        self.assertEqual(seq.sequence, b"MYYK")

    def test_id_property(self):
        seq = Sequence(b"test", b"MYYK")
        self.assertEqual(seq.id, b"test")

    def test_buffer_export(self):
        seq = Sequence(b"test", b"MYYK")
        mem = memoryview(seq)
        self.assertEqual(mem.shape[0], len(seq.sequence))
        self.assertEqual(mem[1], mem[2])
        self.assertTrue(mem.readonly)
        b = mem.tobytes()
        self.assertEqual(len(b), 4)

    def test_error_empty(self):
        with self.assertRaises(ValueError):
            seq = Sequence(b"", b"")


class TestGappedSequence(unittest.TestCase):

    def test_sequence_property(self):
        seq1 = GappedSequence(b"test", b"MY-YK")
        self.assertEqual(seq1.sequence, b"MY-YK")
        seq2 = GappedSequence(b"test", b"--MYYK")
        self.assertEqual(seq2.sequence, b"--MYYK")
