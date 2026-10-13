# Changelog

All notable changes to MSPtools are documented here.

## [Unreleased]

### Fixed

- Preserve physical reference channels as zero log2 ratios during MAD normalization instead of turning them into missing values.

## [1.0.1] - 2026-10-12

### Fixed

- Read plex and channel assignments from OpenMS abundance column headers when the older CSV metadata comment is absent.

## [1.0.0] - 2026-10-11

### Added

- Command line utilities for mzML, idXML, featureXML, consensusXML, protein quantification CSV, and variant peptide annotation.
- Python package distribution, automated tests, and Docker image publishing.

[1.0.0]: https://github.com/TheBoutrosLab/package-MSPtools/releases/tag/v1.0.0
[1.0.1]: https://github.com/TheBoutrosLab/package-MSPtools/compare/v1.0.0...v1.0.1
[unreleased]: https://github.com/TheBoutrosLab/package-MSPtools/compare/v1.0.1...HEAD
