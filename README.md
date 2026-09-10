# Intensities comparator

Compare experimental line intensities from a CSV file against calculated linestrengths from a Fortran `fort.36` output file, and print a comparison table with relative differences.

## Overview

Spectroscopic transition intensities are commonly produced by Fortran programs that write a large line list file (here `fort.36`) containing quantum numbers, transition frequencies, and calculated intensities. This tool matches each *experimental* spectral line to its corresponding *calculated* entry and reports how well the intensities agree.

The matching is driven by quantum-number labels in the experimental CSV (e.g. `P10e`, `R0e`) plus a transition-frequency check.

## Features

- Pure Python 3 standard library — no external dependencies.
- Fast lookup over very large line lists (millions of lines) via an in-memory quantum-number index.
- Branch (`P`/`R`) and symmetry (`e`/`f`) aware matching.
- Hard frequency and intensity tolerances, so transitions that do not genuinely correspond to an experimental line are rejected rather than silently mis-matched.
- Prints a formatted comparison table and summary statistics.
- Tested on N2O calculations.

## Usage

```bash
python compare_intensities.py [options]
```

### Command-line options

| Option          | Default                     | Description |
|-----------------|-----------------------------|-------------|
| `--exp_file`    | `input/Intensities_exp.csv` | Path to the experimental CSV file |
| `--calc_file`   | `input/fort.36`             | Path to the calculated line list (fort.36-like) |
| `--freq_tol`    | `0.01`                      | Absolute frequency tolerance in cm⁻¹ |
| `--int_tol`     | `0.5`                       | Relative intensity tolerance |
| `-h`, `--help`  |                             | Show help and exit |

### Examples

```bash
# Use the defaults
python compare_intensities.py

# Specify both input files
python compare_intensities.py \
    --exp_file input/Intensities_exp_5000.csv \
    --calc_file input/fort.36.D2n

# Loose frequency tolerance
python compare_intensities.py --freq_tol 0.05
```

## Input file formats

### Experimental CSV (`--exp_file`)

Three comma-separated columns: `label, frequency, intensity`.

```
P10e,6364.329281,1.130E-24,
P9e ,6365.286733,1.080E-24,
R0e ,6374.132400,1.394E-25,
```

- **Label** — spectroscopic notation `<branch><J><symmetry>`, for example:
  - `P10e` — a P-branch line starting from rotational level *J* = 10, e-symmetry.
  - `R0e` — an R-branch line starting from *J* = 0, e-symmetry.
  - Branch `P` means *J* decreases by 1 (final *J* = *J* − 1); branch `R` means *J* increases by 1 (final *J* = *J* + 1).
  - Symmetry `e` maps to column 1 = `0` in the line list; `f` maps to `1`.
- **Frequency** — transition frequency in cm⁻¹. Must match the 10th column of the corresponding line list entry.
- **Intensity** — absolute intensity, typically in cm/molecule or similar. Matches the 12th column of the line list.

A trailing empty column (as above) is tolerated.

### Calculated line list (`--calc_file`)

Whitespace-separated Fortran-formatted table with 14 columns. Only the columns below are used:

| Column (1-based) | Meaning                                              |
|------------------|------------------------------------------------------|
| 1                | Symmetry (`0` = e, `1` = f)                          |
| 2                | Final rotational quantum number *J*                  |
| 5                | Initial rotational quantum number *J*                |
| 10               | Transition frequency (cm⁻¹)                          |
| 12               | Calculated intensity                                 |

Example line:

```
  0    1  0  61   2  0   1  6374.133366     2.514086  6371.619280   0.60606389E-07  0.287E-24  0.342E-06  0.164E-02
```

Files like `fort.36` are supported via `--calc_file`.

## Matching algorithm

For each experimental CSV entry, the tool:

1. **Parses the label** into branch (`P`/`R`), initial *J*, and symmetry (`e`/`f`), and derives the final *J* from the branch type.
2. **Builds a candidate set** from the line list: all entries whose (symmetry, final *J*, initial *J*) match — i.e. columns (1, 2, 5).
3. **Hard filter — frequency**: keeps only candidates with `|col10 − freq_exp| ≤ freq_tol` (absolute, cm⁻¹).
4. **Hard filter — intensity**: keeps only candidates with `|col12 − I_exp| / I_exp ≤ int_tol` (relative).
5. **Scores the survivors** and picks the one with the smallest combined score:
   ```
   score = |col10 − freq_exp| / freq_exp + |col12 − I_exp| / I_exp
   ```
6. **Reports the relative intensity difference** as a percentage:
   ```
   Rel Diff (%) = (I_exp / I_calc − 1) × 100
   ```

Entries that fail either hard filter are reported as unmatched with a diagnostic warning, rather than silently matched to a wrong transition.

### Tuning the tolerances

- Genuine frequency offsets between the experimental list and the calculated list for N2O are typically below ~0.006 cm⁻¹, so the default `--freq_tol 0.01` works across line lists spanning 4 → 7200 cm⁻¹.
- If your experimental intensities are in different units or scaling, increase `--int_tol` (e.g. `--int_tol 1.0` for up to 100% relative difference).

## Example output

```
====================================================================================================
COMPARISON TABLE: Experimental vs Calculated Intensities
====================================================================================================
Label      Freq (exp)      Freq (calc)     Freq diff    Int (exp)       Int (calc)      Rel Diff (%)
----------------------------------------------------------------------------------------------------
P10e       6364.329281     6364.330211     -0.000930    1.130E-24       1.140E-24       -0.8772
P9e        6365.286733     6365.287666     -0.000933    1.080E-24       1.070E-24        0.9346
R0e        6374.132400     6374.133339     -0.000939    1.394E-25       1.450E-25       -3.8621
...
====================================================================================================

Summary:
  Successfully matched: 20/20 entries
  Average absolute relative difference: 1.9708%
  Min absolute relative difference: 0.3833%
  Max absolute relative difference: 9.5766%
```

## Requirements

- Python 3.7+ (uses only the standard library).

## Project layout

```
compare_intensities.py      Main script
input/
    Intensities_exp_*.csv   Experimental spectra (various line lists)
    fort.36.D2n, fort.36.D2nB   Calculated line lists
```