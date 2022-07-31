import unittest

from .. import Sequence


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

    def test_buffer_export_empty(self):
        seq = Sequence(b"", b"")
        mem = memoryview(seq)
        self.assertEqual(mem.shape[0], 0)
        self.assertTrue(mem.readonly)
        b = mem.tobytes()
        self.assertEqual(len(b), 0)
