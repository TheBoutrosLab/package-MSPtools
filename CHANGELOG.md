# Changelog

All notable changes to MSPtools are documented here.

## [Unreleased]

## [1.1.0] - 2026-11-08

### Added

- Normalize merged scan and linked multi-plex TMT maps, with optional per-plex median reference-area scaling and a report.

### Changed

- Store linked-output reference intensities as `reference_intensity::<plex_id>`; retain the unsuffixed `reference_intensity` only for features containing one plex. Readers of multi-plex output must use the plex-specific key.

### Fixed

- Reject linked-layout inputs without per-plex MS1 area metadata and prevent scaling reports from overwriting consensusXML files.

## [1.0.3] - 2026-11-07

### Fixed

- Match fixZeroPlex peptide CSV charges to identified peptide hits and report missing mappings clearly.

## [1.0.2] - 2026-10-13

### Added

- Add an optional MADNormalizer quasi-intensity output using each feature's reference intensity.
- Preserve measured zero sample intensities in quasi-intensity output when the reference is valid.
- Print a MAD normalization tally for processed features, finite sample results, and physical reference states.

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
[1.0.2]: https://github.com/TheBoutrosLab/package-MSPtools/compare/v1.0.1...v1.0.2
[1.0.3]: https://github.com/TheBoutrosLab/package-MSPtools/compare/v1.0.2...v1.0.3
[1.1.0]: https://github.com/TheBoutrosLab/package-MSPtools/compare/v1.0.3...v1.1.0
[unreleased]: https://github.com/TheBoutrosLab/package-MSPtools/compare/v1.1.0...HEAD
