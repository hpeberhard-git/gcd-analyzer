# GCD Analyzer - HLA Typing Resolution Quality Analysis

A Perl tool for analyzing HLA (Human Leukocyte Antigen) typing data from ION CSV files to determine the quality of typing resolution. This tool is specifically designed to identify donors with incomplete allele-level resolution or serology-only typing.

## Features

- **Allele-level resolution analysis**: Identifies donors lacking full allele-level HLA typing
- **Serology mode**: Detects donors with serology-only loci (no DNA typing)
- **Age-based statistics**: Breaks down results by donor age
- **Multiple donor types**: Supports both D-Donors and C-Donors
- **Comprehensive reports**: Generates summary CSV files and console output
- **Debug mode**: Detailed information for each donor with incomplete typing

## Requirements

- Perl 5.x with standard modules:
  - `strict`
  - `warnings`
  - `Getopt::Long`

## Installation

```bash
# Clone or download this repository
git clone <repository-url>
cd gcd-analyzer

# Make the script executable (Linux/macOS)
chmod +x analyze_gcd.pl
```

## Usage

### Basic Usage

```bash
# Analyze files in the current directory
./analyze_gcd.pl

# Analyze files in a specific directory
./analyze_gcd.pl --dir /path/to/data

# Use serology mode
./analyze_gcd.pl --dir examples/input --mode serology

# Enable debug output
./analyze_gcd.pl --dir examples/input --debug
```

### Command-line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--dir <directory>` | Input directory containing ION-*-[DC].csv files | Current directory (`.`) |
| `--mode allele` | Require allele-level resolution | `allele` |
| `--mode serology` | Count donors with serology-only loci | - |
| `--debug` | Show detailed debug information for each donor | Off |
| `--help` or `-h` | Display help message | - |

### Example with Test Data

```bash
# Run analysis on the included example files
./analyze_gcd.pl --dir examples/input --mode allele

# View the generated output files
ls examples/input/analysis_results/
```

## Input File Format

The tool expects CSV files with the following naming convention:
- `ION-<number>-D.csv` for D-Donors
- `ION-<number>-C.csv` for C-Donors

### Required CSV Columns

- `ID`: Donor identifier
- `BIRTH_DATE`: Donor birth date in YYYY-MM-DD format
- **DNA typing fields**: `DNA_A1`, `DNA_A2`, `DNA_B1`, `DNA_B2`, `DNA_C1`, `DNA_C2`, `DRB11`, `DRB12`, `DQB11`, `DQB12`
- **Serology fields** (for serology mode): `A1`, `A2`, `B1`, `B2`, `C1`, `C2`, `DR1`, `DR2`, `DQ1`, `DQ2`

## Output

### Console Output

The tool displays:
1. Per-file statistics (by age group)
2. Summary for D-Donors
3. Summary for C-Donors
4. Combined summary (D + C)

Example output:
```
Processing: examples/input/ION-9999-D.csv
Age              Total   No Allele (SER)      Percent
--------------------------------------------------------------------------------
24                   1               0          0.0%
27                   1               1        100.0%
...
```

### CSV Output Files

All output files are saved in `<input_dir>/analysis_results/`:

- `ION-<number>-<type>.csv`: Per-ION analysis results
- `summary_D-Donors.csv`: Summary for all D-Donors
- `summary_C-Donors.csv`: Summary for all C-Donors
- `summary_Overall.csv`: Combined summary

CSV format:
```csv
Age,Total,all_allele_high,Percent
24,1,0,0.0
27,1,1,100.0
...
```

## Analysis Modes

### Allele Mode (Default)

In allele mode, the tool checks if all required HLA loci have allele-level resolution. A valid allele-level typing follows the pattern: `XX:YY:ZZ[:AA][G]` (e.g., `01:01:01G`, `02:01:01`).

**Criteria for "No Allele-Level Resolution":**
- First allele of any required locus is missing or not in allele-level format
- Second allele is present but not in allele-level format (homozygous typing is allowed with one empty allele)

### Serology Mode

In serology mode, the tool identifies donors who have serological typing at one or more loci but lack DNA typing for those same loci.

**Criteria for "Serology-Only":**
- Serological typing present (e.g., `A1` or `A2` has a value)
- Corresponding DNA fields are both empty (e.g., `DNA_A1` and `DNA_A2` are empty)

## Examples

### Example 1: Standard Analysis

```bash
./analyze_gcd.pl --dir examples/input
```

This will analyze all ION CSV files in `examples/input` and check for donors lacking allele-level resolution in any of the required loci.

### Example 2: Serology Analysis with Debug Output

```bash
./analyze_gcd.pl --dir examples/input --mode serology --debug
```

This will:
1. Identify donors with serology-only typing
2. Display detailed information for each affected donor
3. Generate summary statistics

### Example 3: Using the Test Data

The repository includes test data in `examples/input/`:
- `ION-8888-D.csv`: Test file with various typing scenarios
- `ION-9999-D.csv`: Additional test cases

Expected output files are provided in `examples/output/` for comparison.

## Excluded IONs

The tool automatically excludes certain ION codes from analysis. Currently excluded:
- **8405** (KR): Suspected to report serology-only data for all donors despite DNA-typed samples

You can modify the exclusion list by editing the `%exclude_ions` hash in the script.

## Age Calculation

Donor ages are calculated based on the `BIRTH_DATE` field and the current system date. Donors with missing or invalid birth dates are grouped under "unknown" age.

## Troubleshooting

### No CSV files found
- Ensure your input directory contains files matching the pattern `ION-*-[DC].csv`
- Check that you're specifying the correct directory with `--dir`

### Field not found in file
- Verify that your CSV files contain all required columns
- Check that column names match exactly (case-sensitive)

### Permission denied
- On Linux/macOS, ensure the script is executable: `chmod +x analyze_gcd.pl`
- Ensure you have write permissions for the output directory

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT License

## Author

Hans-Peter Eberhard

---

# Deutsch

## Beschreibung

Ein Perl-Tool zur Analyse von HLA-Typisierungsdaten aus ION-CSV-Dateien, um die Qualität der Typisierungsauflösung zu bestimmen. Das Tool identifiziert Spender mit unvollständiger Allel-Level-Auflösung oder ausschließlich serologischer Typisierung.

## Verwendung

```bash
# Grundlegende Verwendung
./analyze_gcd.pl --dir pfad/zu/daten

# Serologischer Modus
./analyze_gcd.pl --dir examples/input --mode serology

# Mit Debug-Ausgabe
./analyze_gcd.pl --dir examples/input --debug

# Hilfe anzeigen
./analyze_gcd.pl --help
```

## Ausgabe

Die Analyse erzeugt:
- Konsolenausgabe mit Statistiken nach Altersgruppen
- CSV-Dateien im Verzeichnis `<eingabe_verzeichnis>/analysis_results/`
- Zusammenfassungen für D-Spender, C-Spender und Gesamt

Weitere Details finden Sie in der englischen Dokumentation oben.
