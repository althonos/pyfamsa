# 🐍🧮 PyFAMSA [![Stars](https://img.shields.io/github/stars/althonos/pyfamsa.svg?style=social&maxAge=3600&label=Star)](https://github.com/althonos/pyfamsa/stargazers)

*[Cython](https://cython.org/) bindings and Python interface to [FAMSA](https://github.com/refresh-bio/FAMSA), an algorithm for ultra-scale multiple sequence alignments.*

[![Actions](https://img.shields.io/github/actions/workflow/status/althonos/pyfamsa/test.yml?branch=main&logo=github&style=flat-square&maxAge=300)](https://github.com/althonos/pyfamsa/actions)
[![Coverage](https://img.shields.io/codecov/c/gh/althonos/pyfamsa?style=flat-square&maxAge=3600&logo=codecov)](https://codecov.io/gh/althonos/pyfamsa/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pyfamsa.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pyfamsa)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pyfamsa?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pyfamsa)
[![AUR](https://img.shields.io/aur/version/python-pyfamsa?logo=archlinux&style=flat-square&maxAge=3600)](https://aur.archlinux.org/packages/python-pyfamsa)
[![Wheel](https://img.shields.io/pypi/wheel/pyfamsa.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pyfamsa/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pyfamsa.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pyfamsa/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pyfamsa.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pyfamsa/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyfamsa/)
[![Mirror](https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400)](https://git.embl.de/larralde/pyfamsa/)
[![Issues](https://img.shields.io/github/issues/althonos/pyfamsa.svg?style=flat-square&maxAge=600)](https://github.com/althonos/pyfamsa/issues)
[![Docs](https://img.shields.io/readthedocs/pyfamsa/latest?style=flat-square&maxAge=600)](https://pyfamsa.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/althonos/pyfamsa/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/badge/dynamic/json?style=flat-square&color=303f9f&maxAge=86400&label=downloads&query=%24.total_downloads&url=https%3A%2F%2Fapi.pepy.tech%2Fapi%2Fprojects%2Fpyfamsa)](https://pepy.tech/project/pyfamsa)

***⚠️ This package is based on FAMSA 2.***

## 🗺️ Overview

[FAMSA](https://github.com/refresh-bio/FAMSA) is a method published in
2016 by Deorowicz *et al.* for large-scale multiple sequence alignments.
It uses state-of-the-art time and memory optimizations as well as a fast
guide tree heuristic to reach very high performance and accuracy.

PyFAMSA is a Python module that provides bindings to [FAMSA](https://github.com/refresh-bio/FAMSA)
using [Cython](https://cython.org/). It implements a user-friendly, Pythonic
interface to align protein sequences using different parameters and access
results directly. It interacts with the FAMSA library interface, which has
the following advantages:

- **single dependency**: pyfamsa is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  FAMSA binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the FAMSA CLI using a
  sub-process and temporary files.
- **friendly interface**: The different guide tree build methods and
  heuristics can be selected from the Python code with a simple keyword
  argument when configuring a new [`Aligner`](https://pyfamsa.readthedocs.io/en/stable/api/aligner.html#pyfamsa.Aligner).

## 🔧 Installing

PyFAMSA can be installed directly from [PyPI](https://pypi.org/project/pyfamsa/),
which hosts some pre-built wheels for the x86-64 architecture (Linux/OSX)
and the Aarch64 architecture (Linux only), as well as the code required to
compile from source with Cython:
```console
$ pip install pyfamsa
```

<!-- Otherwise, pyfamsa is also available as a [Bioconda](https://bioconda.github.io/)
package:
```console
$ conda install -c bioconda pyfamsa
``` -->

Otherwise, have a look at the [Installation page](https://pyfamsa.readthedocs.io/en/stable/install.html) of the [online documentation](https://pyfamsa.readthedocs.io/)

## 💡 Example

Let's create some sequences in memory, align them using the UPGMA method,
(without any heuristic), and simply print the alignment on screen:

```python
from pyfamsa import Aligner, Sequence

sequences = [
    Sequence(b"Sp8",  b"GLGKVIVYGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII"),
    Sequence(b"Sp10", b"DPAVLFVIMLGTITKFSSEWFFAWLGLEINMMVII"),
    Sequence(b"Sp26", b"AAAAAAAAALLTYLGLFLGTDYENFAAAAANAWLGLEINMMAQI"),
    Sequence(b"Sp6",  b"ASGAILTLGIYLFTLCAVISVSWYLAWLGLEINMMAII"),
    Sequence(b"Sp17", b"FAYTAPDLLLIGFLLKTVATFGDTWFQLWQGLDLNKMPVF"),
    Sequence(b"Sp33", b"PTILNIAGLHMETDINFSLAWFQAWGGLEINKQAIL"),
]

aligner = Aligner(guide_tree="upgma")
msa = aligner.align(sequences)

for sequence in msa:
      print(sequence.id.decode().ljust(10), sequence.sequence.decode())
```

This should output the following:
```
Sp10       --------DPAVLFVIMLGTIT-KFS--SEWFFAWLGLEINMMVII
Sp17       ---FAYTAPDLLLIGFLLKTVA-TFG--DTWFQLWQGLDLNKMPVF
Sp26       AAAAAAAAALLTYLGLFLGTDYENFA--AAAANAWLGLEINMMAQI
Sp33       -------PTILNIAGLHMETDI-NFS--LAWFQAWGGLEINKQAIL
Sp6        ------ASGAILTLGIYLFTLCAVIS--VSWYLAWLGLEINMMAII
Sp8        ------GLGKVIVYGIVLGTKSDQFSNWVVWLFPWNGLQIHMMGII
```

## 🧶 Thread-safety

`Aligner` objects are thread-safe, and the `align` method is re-entrant. You
could batch process several alignments in parallel using a
[`ThreadPool`](https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.ThreadPool) with a single
aligner object:
```python
import glob
import multiprocessing.pool
import Bio.SeqIO
from pyfamsa import Aligner, Sequence

families = [
    [ Sequence(r.id.encode(), r.seq.encode()) for r in Bio.SeqIO.parse(file, "fasta") ]
    for file in glob.glob("pyfamsa/tests/data/*.faa")
]

aligner = Aligner()
with multiprocessing.pool.ThreadPool() as pool:
    alignments = pool.map(aligner.align, families)
```

<!-- ## ⏱️ Benchmarks -->

## 🔎 See Also

Done with your protein alignment? You may be interested in trimming it: in that
case, you could use the [`pytrimal`](https://github.com/althonos/pytrimal) Python
package, which wraps [trimAl](http://trimal.cgenomics.org/) 2.0. Or perhaps
you want to build a HMM from the alignment? Then maybe have a look at
[`pyhmmer`](https://github.com/althonos/pyhmmer), a Python package which
wraps [HMMER](http://hmmer.org/).

## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the [GitHub issue tracker](https://github.com/althonos/pyfamsa/issues)
if you need to report or ask something. If you are filing in on a bug,
please include as much information as you can about the issue, and try to
recreate the same bug in a simple, easily reproducible situation.


### 🏗️ Contributing

Contributions are more than welcome! See
[`CONTRIBUTING.md`](https://github.com/althonos/pyfamsa/blob/main/CONTRIBUTING.md)
for more details.


## 📋 Changelog

This project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html)
and provides a [changelog](https://github.com/althonos/pyfamsa/blob/main/CHANGELOG.md)
in the [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) format.


## ⚖️ License

This library is provided under the [GNU General Public License v3.0](https://choosealicense.com/licenses/gpl-3.0/). FAMSA is developed by the
[REFRESH Bioinformatics Group](https://refresh-bio.github.io/) and is
distributed under the terms of the GPLv3 as well. See `vendor/FAMSA/LICENSE`
for more information. In addition, FAMSA vendors several libraries for
compatibility, all of which are redistributed with PyFAMSA under their own
terms: `atomic_wait` (MIT License), `mimalloc` (MIT License), `libdeflate`
(MIT License),  Boost (Boost Software License).

*This project is in no way not affiliated, sponsored, or otherwise endorsed
by the [FAMSA authors](https://github.com/refresh-bio). It was developed
by [Martin Larralde](https://github.com/althonos/) during his PhD project
at the [European Molecular Biology Laboratory](https://www.embl.de/) in
the [Zeller team](https://github.com/zellerlab).*
