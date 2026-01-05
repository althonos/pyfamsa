Interfacing with other libraries
================================

PyFAMSA is flexible over its inputs and outputs, and it is easy
to integrate with other Python libraries from the bioinformatics 
ecosystem. This page lists some recipes on how to convert datatypes
between these libraries to facilitate building a larger workflow.

.. note::

    Feel free to open a pull request in the 
    `GitHub repository <https://github.com/althonos/pyfamsa>`_ if you 
    have some recipes you would like to add to the list!


PyHMMER
-------

Converting a `pyhmmer.easel.TextSequence` to a `pyfamsa.Sequence`
to pass to a `pyfamsa.Aligner`:

.. code:: python

    import pyfamsa
    import pyhmmer

    with pyhmmer.easel.SequenceFile("sequences.fasta") as f:
        text_seq = f.read()

    sequence = pyfamsa.Sequence(id=text_seq.id, sequence=text_seq.sequence)  


Converting a `pyfamsa.Alignment` to a `pyhmmer.easel.TextMSA` that
can then be digitized and used with a `pyhmmer.plan7.Builder` to build
a HMM:

.. code:: python

    import pyfamsa
    import pyhmmer

    aligner = pyfamsa.Aligner()
    pyfamsa_alignment = aligner.align(sequences)

    text_msa = pyhmmer.easel.TextMSA(
        name=b"alignment",
        sequences=[
            pyhmmer.easel.TextSequence(name=row.id, sequence=row.sequence)
            for row in pyfamsa_alignment
        ]
    )

Converting a `pyhmmer.easel.TextMSA` loaded from a file or built
from the PyHMMER API into a `pyfamsa.Alignment`, for instance to
perform a profile-profile alignment with `pyfamsa.Aligner.align_profiles`:

.. code:: python

    import pyfamsa
    import pyhmmer

    with pyhmmer.easel.MSAFile("sequences.sto") as f:
        text_msa = f.read()

    pyfamsa_alignment = pyfamsa.Alignment(
        pyfamsa.GappedSequence(id=id_, sequence=sequence)
        for id_, sequence in zip(text_msa.name, text_msa.aligned)
    )


Pytrimal
--------

Converting a `pyfamsa.Alignment` to a `pytrimal.Alignment` for trimming
the generated alignment:

.. code:: python

    import pyfamsa
    import pytrimal

    aligner = pyfamsa.Aligner()
    pyfamsa_alignment = aligner.align(sequences)

    pytrimal_alignment = pytrimal.Alignment(
        names=[row.id for row in pyfamsa_alignment],
        sequences=[row.sequence for row in pyfamsa_alignment]    
    )

    trimmer = pytrimal.AutomaticTrimmer()
    trimmed_alignment = trimmer.trim(pytrimal_alignment)

        
PyCoMSA
-------

Converting a `pyfamsa.Alignment` to `pycomsa.MSA` to write an alignment
computed with FAMSA into a CoMSA file with compression:

.. code:: python

    import pycomsa
    import pyfamsa

    aligner = pyfamsa.Aligner()
    pyfamsa_alignment = aligner.align(sequences)

    pycomsa_msa = pycomsa.MSA(
        id="alignment",
        names=[row.id.decode() for row in pyfamsa_alignment],
        sequences=[row.sequence.decode() for row in pyfamsa_alignment]  
    )

    with pycomsa.open("sequences.msac", "w", format="stockholm") as f:
        f.write(pycomsa_msa)


Biopython
---------

Converting `~Bio.SeqRecord.SeqRecord` objects to `pyfamsa.Sequence`
so they can be aligned:

.. code:: python

    import Bio.SeqIO
    import pyfamsa

    sequences = [
        pyfamsa.Sequence(id=record.id.encode(), sequence=bytes(record.seq))
        for record in Bio.SeqIO.parse("sequences.fa", "fasta"):
    ]
       
    aligner = pyfamsa.Aligner()
    pyfamsa_alignment = aligner.align(sequences)

