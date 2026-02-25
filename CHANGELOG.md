# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-02-25

### Added
- Initial release of GCD Analyzer
- Allele-level resolution analysis mode
- Serology mode for detecting serology-only loci
- Age-based statistics and reporting
- Support for D-Donors and C-Donors
- Debug mode for detailed donor information
- CSV output for all analysis results
- Summary reports (D-Donors, C-Donors, Overall)
- Configurable input directory via `--dir` option
- Comprehensive README with usage examples
- Example input and output files
- MIT License

### Changed
- Default input directory changed from hard-coded path to current directory (`.`)
- Improved help message with examples

### Features
- Automatic exclusion of specific ION codes (currently ION-8405)
- Age calculation based on birth date
- Validation of required CSV columns
- Error handling for missing or malformed data

## [Unreleased]

### Planned
- Configuration file support for excluded IONs
- Additional output formats (JSON, Excel)
- Web-based interface
- Batch processing improvements
- Statistical summaries and visualizations
