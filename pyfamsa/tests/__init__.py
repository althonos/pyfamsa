# noqa: D104

from . import (
    test_aligner,
    test_alignment,
    test_doctest,
    test_sequence,
)


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_aligner))
    suite.addTests(loader.loadTestsFromModule(test_alignment))
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_sequence))
    return suite
