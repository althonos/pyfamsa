Utilities
=========

.. currentmodule:: pyfamsa

Constants
---------

.. data:: MIQS
    :type: scoring_matrices.ScoringMatrix

    The MIQS scoring matrix proposed by Yamada & Tomii (2014), used by default
    in FAMSA for scoring alignments until FAMSA ``v2.3.0`` (PyFAMSA ``v0.6.0``).

.. data:: PFASUM31
    :type: scoring_matrices.ScoringMatrix

    The PFASUM31 scoring matrix proposed by Keul & Hess (2017).

.. data:: PFASUM43
    :type: scoring_matrices.ScoringMatrix

    The PFASUM43 scoring matrix proposed by Keul & Hess (2017) and used by
    default in FAMSA for scoring alignments since FAMSA ``v2.3.0``
    (PyFAMSA ``v0.6.0``).

.. data:: PFASUM60
    :type: scoring_matrices.ScoringMatrix

    The PFASUM60 scoring matrix proposed by Keul & Hess (2017).

.. data:: FAMSA_ALPHABET
    :type: str

    The alphabet used by FAMSA to encode input sequences.


Functions
---------

.. autofunction:: pyfamsa.famsa_info
