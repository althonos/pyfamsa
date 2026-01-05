# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyfamsa/compare/v0.6.2...HEAD


## [v0.6.2] - 2026-01-05
[v0.6.2]: https://github.com/althonos/pyfamsa/compare/v0.6.1...v0.6.2

### Fixed
- `std::fill` and `std::vector` C++ algorithm functions changing signatures in latest Cython version.
- Patch segmentation fault in profile-profile alignments ([#8](https://github.com/althonos/pyfamsa/issues/8), [refresh-bio/FAMSA#61](https://github.com/refresh-bio/FAMSA/issues/61)).
- Ensure build in debug mode on ReadTheDocs. 


## [v0.6.1] - 2025-10-02
[v0.6.1]: https://github.com/althonos/pyfamsa/compare/v0.6.0...v0.6.1

### Added
- Sequence data validation against `Aligner.scoring_matrix.alphabet`.
- Compilation of Limited API wheels for Python 3.11 and later.

### Changed
- Support scoring matrices with alphabets that are subsets of `FAMSA_ALPHABET`.
- Avoid unneeded copy in `Sequence.__init__` and `GappedSequence.__init__`.

### Fixed
- Incorrect handling of extra symbols in `Sequence` and `Aligner`.


## [v0.6.0] - 2025-08-13
[v0.6.0]: https://github.com/althonos/pyfamsa/compare/v0.5.3-post1...v0.6.0

### Added
- `PFASUM31`, `PFASUM43` and `PFASUM60` scoring matrices as `pyfamsa` constants.

### Changed
- Update FAMSA to `v2.4.1`.
- Change default scoring matrix from *MIQS* to *PFASUM41*.
- Drop support for Python 3.7 due to outdated `manylinux` setup.`


## [v0.5.3-post1] - 2025-03-04
[v0.5.3-post1]: https://github.com/althonos/pyfamsa/compare/v0.5.3...v0.5.3-post1

### Fixed
- Extra key in `pyproject.toml` causing build issues with version `0.11.0` of `scikit-build-core`.


## [v0.5.3] - 2024-11-05
[v0.5.3]: https://github.com/althonos/pyfamsa/compare/v0.5.2...v0.5.3

### Added
- Support for Python 3.13.

### Changed
- Use `scikit-build-core` to build package.

### Removed
- Support for Python 3.6.



## [v0.5.2] - 2024-09-19
[v0.5.2]: https://github.com/althonos/pyfamsa/compare/v0.5.1...v0.5.2

### Changed
- Update to FAMSA [`v2.2.3`](https://github.com/refresh-bio/FAMSA/releases/tag/v2.2.3).


## [v0.5.1] - 2024-08-28
[v0.5.1]: https://github.com/althonos/pyfamsa/compare/v0.5.0...v0.5.1

### Fixed
- Unit tests failing on missing on missing optional `importlib-resources` dependency.


## [v0.5.0] - 2024-08-28
[v0.5.0]: https://github.com/althonos/pyfamsa/compare/v0.4.0...v0.5.0

### Added
- Constructor to `GappedSequence` class, taking an identifier and a sequence as `bytes` objects.
- Constructor to `Alignment` class, taking an iterable of `GappedSequence` object.
- `Alignment.copy` implementation.
- Slicing implementation to `Alignment`.
- `Aligner.align_profiles` function to align two profiles ([#5](https://github.com/althonos/pyfamsa/issues/5)).

### Fixed
- Pin supported versions of `scoring-matrices` package to `~=0.2.0`.

### Changed
- Use C++ `shared_ptr` in `GappedSequence` and `Alignment` to avoid copying data when possible.
- Migrate documentation to `pydata-sphinx-theme`.


## [v0.4.0] - 2024-05-06
[v0.4.0]: https://github.com/althonos/pyfamsa/compare/v0.3.2...v0.4.0

### Added
- `scoring-matrices` dependency to handle alternative scoring matrices.
- `scoring_matrix` argument to `Aligner` constructor to use a non-default matrix ([#3](https://github.com/althonos/pyfamsa/issues/3)).

### Fixed
- Use of outdated `importlib.resources` interface in `pyfamsa.tests` package.
- Missing defines for compilation of NEON code on non-Aarch64 Arm platforms.


## [v0.3.2] - 2024-01-27
[v0.3.2]: https://github.com/althonos/pyfamsa/compare/v0.3.1...v0.3.2

### Added
- `pickle` protocol support for `Sequence` objects.

### Fixed
- Disable creation of empty `Sequence` objects to prevent segmentation faults in FAMSA ([#2](https://github.com/althonos/pyfamsa/issues/1)).


## [v0.3.1] - 2023-01-14
[v0.3.1]: https://github.com/althonos/pyfamsa/compare/v0.3.0...v0.3.1

### Fixed
- Disable use of memory-monotonic allocations to fix multithreading errors ([#1](https://github.com/althonos/pyfamsa/issues/1)).


## [v0.3.0] - 2023-07-21
[v0.3.0]: https://github.com/althonos/pyfamsa/compare/v0.2.0...v0.3.0

### Changed
- Bumped Cython dependency to `v3.0`.

### Fixed
- PyPy builds failing on missing `PyInterpreterState_GetID`.


## [v0.2.0] - 2022-11-22
[v0.2.0]: https://github.com/althonos/pyfamsa/compare/v0.1.1...v0.2.0

### Added
- `pyfamsa.famsa_info` function to get version information about the embedded FAMSA version.
- Explicit support for Python 3.11.
- Wheel distributions for MacOS Aarch64 platforms.

### Changed
- Bumped vendored FAMSA to `v2.2.2`.

### Fixed
- `Aligner.build_tree` and `Aligner.align` now accept inputs containing less than two sequences.



## [v0.1.1] - 2022-08-06
[v0.1.1]: https://github.com/althonos/pyfamsa/compare/v0.1.0...v0.1.1

### Added
- MyPy Type stubs for the `pyfamsa._famsa` Cython extension.
- Documentation for the `Aligner.build_tree` method.

### Fixed
- Missing header files for compilation on older platforms.
- Missing define macros for Windows build target.


## [v0.1.0] - 2022-08-05
[v0.1.0]: https://github.com/althonos/pyfamsa/compare/5dca3122...v0.1.0

Initial release.
