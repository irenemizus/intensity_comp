# Intensities comparator

Compare experimental line intensities from a CSV file against calculated linestrengths from a Fortran `fort.36` output file, and print a comparison table with relative differences.

## Overview

Spectroscopic transition intensities are commonly produced by Fortran programs that write a large line list file (here `fort.36`) containing quantum numbers, transition frequencies, and calculated intensities. This tool matches each *experimental* spectral line to its corresponding *calculated* entry and reports how well the intensities agree.

The matching is driven by quantum-number labels in the experimental CSV (e.g. `P10e`, `R0e`) plus a transition-frequency check.

## Features

- Pure Python 3 standard library for the core comparison — no external dependencies required for the default workflow.
- Fast lookup over very large line lists (millions of lines) via an in-memory quantum-number index.
- Branch (`P`/`R`) and symmetry (`e`/`f`) aware matching.
- Hard frequency and intensity tolerances, so transitions that do not genuinely correspond to an experimental line are rejected rather than silently mis-matched.
- Prints a formatted comparison table and summary statistics.
- **Optional per-transition precision labelling** (`--rel_weights`): each transition is labelled by the precision of its initial and final states (`HP (HP/HP)`, `LP (LP/LP)`, `Mixed (HP/LP)`, `Mixed (LP/HP)`), with per-label counts and per-label sd values added to the summary. This mode uses the read-only [`sdcalcs`](#the-sdcalcs-dependency) package and an experimental energy-level file.
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
| `--rel_weights` | *(off)*                     | Switch on the per-transition precision labelling (see [Precision labelling](#precision-labelling--rel_weights)) |
| `--eps`         | `0.05`                      | State-energy tolerance (cm⁻¹) used to match a calculated state to an experimental level; only used with `--rel_weights` |
| `--energy_file` | `input/EnergyLevels_100_n2o.txt` | Path to the experimental energy-level file with uncertainties; only used with `--rel_weights` |
| `--out_file`    | *(none)*                      | If given, also write the comparison table to this file as **tab-separated values** (with a header row) — a format that pastes directly into Excel/other spreadsheets |
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

 # Enable per-transition precision labelling (uses the energy-level file
 # and the sdcalcs package)
 python compare_intensities.py \
     --exp_file input/Intensities_exp_5000.csv \
     --calc_file input/fort.36.D2n \
     --rel_weights \
     --energy_file input/EnergyLevels_100_n2o.txt \
     --eps 0.05

 # Also write the comparison table to a tab-separated file (Excel-friendly);
 # the precision column is included when --rel_weights is on
 python compare_intensities.py \
     --exp_file input/Intensities_exp_5000.csv \
     --calc_file input/fort.36.D2n \
     --rel_weights \
     --out_file table.tsv
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

Whitespace-separated Fortran-formatted table with 14 columns. The columns below are used:

| Column (1-based) | Meaning                                              | Used by |
|------------------|------------------------------------------------------|---------|
| 1                | Constant `0` in the current files; kept as a matching-key prefix | matching |
| 2                | Final rotational quantum number *J*                  | matching |
| 3                | Symmetry of the final state (`0` = e, `1` = f)       | precision labelling |
| 5                | Initial rotational quantum number *J*                | matching |
| 6                | Symmetry of the initial state (`0` = e, `1` = f)     | precision labelling |
| 8                | Final-state energy *E* (cm⁻¹)                        | precision labelling |
| 9                | Initial-state energy *E* (cm⁻¹)                      | precision labelling |
| 10               | Transition frequency (cm⁻¹)                          | matching |
| 12               | Calculated intensity                                 | matching |

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

## Precision labelling (`--rel_weights`)

Each transition connects an **initial** state and a **final** state. The precision of a state is taken from the experimental energy-level file (`--energy_file`) and depends on how many transitions constrain that level: a level is **high-precision (HP)** when its per-transition experimental error `dE = unc · √(#lines) ≤ 0.001 cm⁻¹` and **low-precision (LP)** otherwise. This HP/LP split is computed by the [`sdcalcs`](#the-sdcalcs-dependency) package (`assign_precision_weights`).

When `--rel_weights` is on, every successfully matched transition is additionally labelled by the precision pair of its initial/final states:

| Label              | Meaning                          |
|--------------------|----------------------------------|
| `HP (HP/HP)`       | initial HP, final HP             |
| `LP (LP/LP)`       | initial LP, final LP             |
| `Mixed (HP/LP)`    | initial HP, final LP             |
| `Mixed (LP/HP)`    | initial LP, final HP             |
| `Unresolved (X/Y)` | at least one state could not be matched; `X`/`Y` are the initial/final precision — `HP`, `LP` if that state was found, or `NA` if it was not (e.g. `Unresolved (HP/NA)` = initial HP, final not found) |

The label is printed in a new **`Precision type (Initial/Final)`** column of the comparison table, and the summary gains a *Precision label breakdown* block reporting, for each label, the number of transitions `N` and the (unweighted) standard deviation of the relative intensity differences `sd = √(Σ d² / N)`.

### How a transition is resolved to its two levels

For a matched `fort.36` row the tool takes the final state `(J_final, sym_final, E_final)` (columns 2, 3, 8) and the initial state `(J_initial, sym_initial, E_initial)` (columns 5, 6, 9), and looks each of them up in the energy-level file:

- the **final** level is the experimental level with `J = J_final`, `sym = sym_final` whose energy `E` satisfies `|E − E_final| ≤ --eps`;
- the **initial** level is the experimental level with `J = J_initial`, `sym = sym_initial` whose energy `E` satisfies `|E − E_initial| ≤ --eps`.

When both are found the transition is labelled; when either is missing the transition is kept in the table but marked **`Unresolved (X/Y)`** (where `X`/`Y` are the precision of the initial/final state that *was* found — `HP`/`LP` — or `NA` if that state was not found), and a warning listing the offending state (its `J`, symmetry and calculated energy) is printed after the summary. In the summary breakdown all `Unresolved (…)` variants are grouped under a single **`Unresolved`** total. `--eps` (default `0.05 cm⁻¹`) controls how close a calculated state energy must be to an experimental level.

### Energy-level file format (`--energy_file`)

Whitespace-separated N2O layout, one level per line:

| Column (1-based) | Meaning |
|------------------|---------|
| 1                | vibrational quantum number *v1* |
| 2                | vibrational angular-momentum quantum number *l* |
| 3                | vibrational quantum number *v2* |
| 4                | symmetry (`e`/`f`) |
| 5                | rotational quantum number *J* |
| 6                | energy *E* (cm⁻¹) |
| 7                | experimental uncertainty `unc` (cm⁻¹) |
| 8                | number of constraining transitions `#lines` |

Example:

```
0 0 1 e 0 	0.000000000	0.000e+00	139
0 0 1 e 1 	0.838021281	2.068e-08	335
```

## The sdcalcs dependency

The precision-labelling mode is the only part that depends on the external `sd_calcs` project. Its `sdcalcs` package is imported **read-only** from the sibling directory `../sd_calcs` (relative to this script); nothing in that project is modified. If that directory is not present (or not a sibling), `--rel_weights` exits with an explanatory error, while the default workflow (without `--rel_weights`) needs no external package at all.

## Example output

### Default (without `--rel_weights`)

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

### With `--rel_weights`

The comparison table gains a trailing `Precision type (Initial/Final)` column and the summary gains a per-label breakdown:

```
...
Label      Freq (exp)      Freq (calc)     Freq diff    Int (exp)       Int (calc)      Rel Diff (%)    Precision type (Initial/Final)
------------------------------------------------------------------------------------------------------------------------------------
P10e       6364.329281     6364.330211     -0.000930    1.130E-24       1.140E-24       -0.8772         HP (HP/HP)
P9e        6365.286733     6365.287666     -0.000933    1.080E-24       1.070E-24        0.9346         HP (HP/HP)
R0e        6374.132400     6374.133339     -0.000939    1.394E-25       1.450E-25       -3.8621         Mixed (HP/LP)
P2e        6371.618341     6371.619280     -0.000939    2.859E-25       2.870E-25       -0.3833         Unresolved (HP/NA)
...
====================================================================================================

Summary:
  Successfully matched: 20/20 entries
  Average absolute relative difference: 1.9708%
  Standard deviation: 2.4051%
  Min absolute relative difference: 0.3833%
  Max absolute relative difference: 9.5766%

Precision label breakdown (Initial/Final):
  HP (HP/HP)          : N = 8, sd = 1.7010% (over 8 with a valid relative difference)
  LP (LP/LP)          : N = 0
  Mixed (HP/LP)       : N = 11, sd = 3.4441% (over 11 with a valid relative difference)
  Mixed (LP/HP)       : N = 0
  Unresolved          : N = 1, sd = 0.3833% (over 1 with a valid relative difference)

WARNING: 1 transition(s) could not be resolved to a pair of experimental levels
(no level within 0.050 cm^-1 of the calculated state energy):
  P2e (freq=6371.6183 cm^-1): no experimental final state (J=1, sym=e, E_calc=6374.1333 cm^-1)
```

## Writing the table to a file (`--out_file`)

Pass `--out_file <path>` to also write the comparison table to a file as **tab-separated values** (TSV): a header row followed by one tab-separated row per transition, mirroring the on-screen table. TSV is chosen because it pastes cleanly into Excel/Google Sheets/Numbers (the tab is auto-detected as the column separator).

- Columns: `Label`, `Freq (exp)`, `Freq (calc)`, `Freq diff`, `Int (exp)`, `Int (calc)`, `Rel Diff (%)` — plus `Precision type (Initial/Final)` when `--rel_weights` is on.
- Intensities are written in scientific notation (e.g. `1.130E-24`); frequencies and the relative difference use fixed-point decimals.
- Transitions with no match (and an `inf` relative difference) are written with `N/A` in the affected cells.

Example (`--rel_weights --out_file table.tsv`):

```
Label	Freq (exp)	Freq (calc)	Freq diff	Int (exp)	Int (calc)	Rel Diff (%)	Precision type (Initial/Final)
P10e	6364.329281	6364.330211	-0.000930	1.130E-24	1.140E-24	-0.8772	HP/HP
P9e	6365.286733	6365.287666	-0.000933	1.080E-24	1.070E-24	0.9346	HP/HP
R0e	6374.132400	6374.133339	-0.000939	1.394E-25	1.450E-25	-3.8621	HP/LP
P2e	6371.618341	6371.619280	-0.000939	2.859E-25	2.870E-25	-0.3833	HP/NA
```

Only the comparison table is written (not the summary/warnings), so the file is ready to paste into a spreadsheet.

## Requirements

- Python 3.7+ (the default workflow uses only the standard library).
- For the `--rel_weights` precision-labelling mode: the read-only `sdcalcs` package from the sibling `../sd_calcs` project (imported automatically; no installation needed).

## Project layout

```
compare_intensities.py      Main script
input/
    Intensities_exp_*.csv   Experimental spectra (various line lists)
    fort.36.D2n, fort.36.D2nB   Calculated line lists
    EnergyLevels_100_n2o.txt    Experimental energy levels + uncertainties (--rel_weights)
../sd_calcs/                Read-only dependency for --rel_weights (sdcalcs package)
```