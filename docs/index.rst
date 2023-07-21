PyFAMSA |Stars|
===============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyfamsa.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyfamsa/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `FAMSA <https://github.com/refresh-bio/FAMSA>`_,
*an algorithm for ultra-scale multiple sequence alignments.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyfamsa/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyfamsa/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyfamsa/branch/main.svg?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyfamsa/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyfamsa.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyfamsa

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyfamsa?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyfamsa

.. |AUR| image:: https://img.shields.io/aur/version/python-pyfamsa?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyfamsa

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyfamsa?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyfamsa/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyfamsa.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyfamsa/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyfamsa.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyfamsa/#files

.. |License| image:: https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400
   :target: https://choosealicense.com/licenses/gpl-3.0/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyfamsa/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pyfamsa/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyfamsa.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyfamsa/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyfamsa?style=flat-square&maxAge=3600
   :target: http://pyfamsa.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyfamsa/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyfamsa?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyfamsa



Overview
--------

`FAMSA <https://github.com/refresh-bio/FAMSA>`_ is a method published in
2016 by Deorowicz *et al.* for large-scale multiple sequence alignments.
It uses state-of-the-art time and memory optimizations as well as a fast
guide tree heuristic to reach very high performance and accuracy.

``pyfamsa`` is a Python module that provides bindings to FAMSA using the
`Cython <https://cython.org/>`_ language. It implements a user-friendly,
Pythonic interface to align protein sequences using different parameters and
access results directly. It interacts with the FAMSA library interface, which
has the following advantages:

- **single dependency**: ``pyfamsa`` is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  FAMSA binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the FAMSA CLI using a
  sub-process and temporary files.
- **friendly interface**: The different guide tree build methods and
  heuristics can be selected from the Python code with a simple keyword
  argument when configuring a new `Aligner`.


Setup
-----

Run ``pip install pyfamsa`` in a shell to download the latest release and all
its dependencies from PyPi, or have a look at the
:doc:`Installation page <install>` to find other ways to install ``pyfamsa``.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Contributing <contributing>
   API Reference <api/index>
   Changelog <changes>


License
-------

This library is provided under the `GNU General Public License v3.0 <https://choosealicense.com/licenses/gpl-3.0/>`_.

FAMSA is developed by the `REFRESH Bioinformatics Group <https://refresh-bio.github.io/>`_
and is distributed under the terms of the GPLv3 as well.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original FAMSA authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
