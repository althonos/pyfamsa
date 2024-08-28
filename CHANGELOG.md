# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyfamsa/compare/v0.5.1...HEAD


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
