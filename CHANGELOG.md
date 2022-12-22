# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased]
[Unreleased]: https://github.com/althonos/pyfamsa/compare/v0.2.0...HEAD


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
