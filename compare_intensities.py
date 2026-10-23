#!/usr/bin/env python3
"""
Compare experimental intensities from CSV with calculated values from fort.36.

Usage:
    python compare_intensities.py [options]

Options:
    --exp_file  Path to the experimental CSV file (default: input/Intensities_exp.csv)
    --calc_file Path to the fort.36-like calculated file (default: input/fort.36)
    --freq_tol  Absolute frequency tolerance in cm^-1 (default: 0.01)
    --int_tol   Relative intensity tolerance (default: 0.5)

    --rel_weights  Switch on the per-transition precision labelling. Each transition
                   is labelled by the precision of its initial and final states, where a
                   state is high-precision (HP) when its per-transition experimental error
                   dE = unc * sqrt(#lines) <= 0.001 cm^-1 and low-precision (LP) otherwise.
                   The precision of the states is taken from the experimental energy-level
                   file (see --energy_file) via the sdcalcs API.
    --eps         Tolerance (cm^-1) used to match a calculated state energy (from fort.36)
                   to an experimental level; only used with --rel_weights (default: 0.05)
    --energy_file Path to the experimental energy-level file with uncertainties, in the
                   N2O layout 'v1 l v2 sym J E unc #lines'
                   (default: input/EnergyLevels_100_n2o.txt); only used with --rel_weights
"""
import bisect
import math
import os
import re
import sys
import argparse
from collections import defaultdict


def _import_sdcalcs():
    """
    Import the read-only 'sdcalcs' package from the sibling 'sd_calcs' project.

    The two projects are expected to be siblings (../sd_calcs relative to this file).
    The path is added to sys.path on the fly so the package can be imported without
    modifying or installing the sd_calcs project.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidate = os.path.normpath(os.path.join(script_dir, '..', 'sd_calcs'))
    if os.path.isdir(os.path.join(candidate, 'sdcalcs')) and candidate not in sys.path:
        sys.path.insert(0, candidate)
    import sdcalcs
    return sdcalcs


def build_level_index(energy_file, dE_th=0.001):
    """
    Parse the experimental energy-level file and classify every level as
    high-precision (HP) or low-precision (LP) using the sdcalcs API.

    Returns a dict mapping (J, sym_char) -> list of (E, is_hp) sorted by E, where
    sym_char is 'e' or 'f'.
    """
    sdcalcs = _import_sdcalcs()
    # The level file only needs the levels actually present; a generous J range
    # guarantees nothing is filtered out (filter_Jlist_lines keeps J in J_list).
    J_list = [str(J) for J in range(0, 200)]
    states_exp = sdcalcs.ExpFormatN2O(energy_file, J_list, 4).parse_file()
    sdcalcs.assign_precision_weights(states_exp, dE_th, 1.0)

    lev_index = defaultdict(list)
    for st in states_exp:
        # ExpFormatN2O maps 'e' -> SymType.A1 (1) and 'f' -> SymType.A2 (2)
        sym_char = 'e' if st.sym == 1 else 'f'
        lev_index[(st.J, sym_char)].append((st.E, st.is_hp))
    for key in lev_index:
        lev_index[key].sort(key=lambda x: x[0])
    return lev_index


def find_level(lev_index, J, sym_char, E_target, eps):
    """
    Find the experimental level with the given (J, symmetry) whose energy is closest
    to E_target. Returns (abs_diff, is_hp, E) when the closest level lies within
    `eps` of E_target, otherwise None.
    """
    candidates = lev_index.get((J, sym_char), [])
    if not candidates:
        return None
    energies = [c[0] for c in candidates]
    i = bisect.bisect_left(energies, E_target)
    best = None
    for j in (i - 1, i):
        if 0 <= j < len(candidates):
            E, is_hp = candidates[j]
            d = abs(E - E_target)
            if best is None or d < best[0]:
                best = (d, is_hp, E)
    if best is None or best[0] > eps:
        return None
    return best


def _match_precision(cols, lev_index, eps):
    """
    Given a matched fort.36 row (list of string columns), resolve its initial and
    final states to experimental levels and return the precision label.

    fort.36 column layout (1-indexed as documented):
        2 -> J_final, 3 -> sym_final (0=e, 1=f), 5 -> J_initial, 6 -> sym_initial,
        8 -> E_final, 9 -> E_initial, 10 -> frequency.

    Returns (label, details):
        label is one of 'HP (HP/HP)', 'LP (LP/LP)', 'Mixed (HP/LP)', 'Mixed (LP/HP)'
        when both states are resolved, or 'Unresolved (<init>/<fin>)' when at least one
        state cannot be matched to an experimental level within `eps` — where each of
        <init>/<fin> is 'HP', 'LP' (when that state was found) or 'NA' (when it was not
        found), e.g. 'Unresolved (HP/NA)'.
        details is a dict describing the calculated state energies (used for warnings).
    """
    def _state_str(found):
        if found is None:
            return 'NA'
        return 'HP' if found[1] else 'LP'

    try:
        J_final = int(cols[1]); sym_final = 'e' if int(cols[2]) == 0 else 'f'
        J_initial = int(cols[4]); sym_initial = 'e' if int(cols[5]) == 0 else 'f'
        E_final = float(cols[7]); E_initial = float(cols[8])
    except (ValueError, IndexError):
        return 'Unresolved (NA/NA)', {'error': 'cannot parse fort.36 row columns',
                                      'fin': None, 'ini': None}

    fin = find_level(lev_index, J_final, sym_final, E_final, eps)
    ini = find_level(lev_index, J_initial, sym_initial, E_initial, eps)

    details = {
        'J_final': J_final, 'sym_final': sym_final, 'E_final': E_final, 'fin': fin,
        'J_initial': J_initial, 'sym_initial': sym_initial, 'E_initial': E_initial, 'ini': ini,
    }
    init_str = _state_str(ini)
    fin_str = _state_str(fin)
    if fin is None or ini is None:
        return f'Unresolved ({init_str}/{fin_str})', details

    if init_str == 'HP' and fin_str == 'HP':
        label = 'HP (HP/HP)'
    elif init_str == 'LP' and fin_str == 'LP':
        label = 'LP (LP/LP)'
    elif init_str == 'HP' and fin_str == 'LP':
        label = 'Mixed (HP/LP)'
    else:
        label = 'Mixed (LP/HP)'
    return label, details


def _print_label_stats(name, group):
    """Print the count and unweighted sd of the relative intensity differences for a
    group of comparison rows (row[5] = rel_diff)."""
    n = len(group)
    if n == 0:
        print(f"  {name:<20}: N = 0")
        return
    valid = [r[5] for r in group if r[5] is not None and r[5] != float('inf')]
    if valid:
        sd_g = math.sqrt(sum(d ** 2 for d in valid) / len(valid))
        print(f"  {name:<20}: N = {n}, sd = {sd_g:.4f}% "
              f"(over {len(valid)} with a valid relative difference)")
    else:
        print(f"  {name:<20}: N = {n}")


def parse_label(label):
    """
    Parse a CSV label like 'P10e' or 'R0e' into components.
    
    Returns: (branch, J_initial, symmetry)
    - branch: 'P' or 'R'
    - J_initial: int, initial rotational quantum number
    - symmetry: 'e' or 'f'
    """
    label = label.strip()
    match = re.match(r'^([PR])(\d+)(e|f)$', label, re.IGNORECASE)
    if not match:
        raise ValueError(f"Cannot parse label: '{label}'")
    
    branch = match.group(1).upper()
    J_initial = int(match.group(2))
    symmetry = match.group(3).lower()
    
    return branch, J_initial, symmetry


def calculate_J_final(branch, J_initial):
    """
    Calculate final J quantum number based on branch type.
    
    P branch: J decreases by 1 (J_final = J_initial - 1)
    R branch: J increases by 1 (J_final = J_initial + 1)
    """
    if branch == 'P':
        return J_initial - 1
    elif branch == 'R':
        return J_initial + 1
    else:
        raise ValueError(f"Unknown branch type: {branch}")


def parse_csv(csv_file):
    """
    Parse the experimental intensities CSV file.
    
    Returns: list of (label, frequency, intensity)
    """
    data = []
    with open(csv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Split by comma
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 3:
                continue
            
            label = parts[0].strip()
            try:
                frequency = float(parts[1])
                intensity = float(parts[2])
                data.append((label, frequency, intensity))
            except ValueError:
                continue
    
    return data


def parse_fort36(fort36_file):
    """
    Parse the fort.36 file and build an index for fast lookup.
    
    Returns: dict mapping (col1, col2, col5) -> list of (col10, col12, full_cols)
    """
    index = defaultdict(list)
    
    with open(fort36_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            cols = line.split()
            if len(cols) < 14:
                continue
            
            try:
                col1 = int(cols[0])  # symmetry: 0=e, 1=f
                col2 = int(cols[1])  # J_final
                col5 = int(cols[4])  # J_initial
                col10 = float(cols[9])  # transition frequency
                col12 = float(cols[11])  # intensity
            except (ValueError, IndexError):
                continue
            
            key = (col1, col2, col5)
            index[key].append((col10, col12, cols))
    
    return index


def relative_intensity_diff(calc_int, exp_int):
    """
    Relative intensity difference |calc - exp| / exp.
    
    Returns float('inf') when the experimental intensity is zero but the
    calculated one is not (a perfect mismatch), and 0.0 when both are zero.
    """
    if exp_int == 0:
        return 0.0 if calc_int == 0 else float('inf')
    return abs(calc_int - exp_int) / exp_int


def find_matching_entry(J_initial, J_final, symmetry, frequency, intensity, 
                        fort36_index, freq_tol, int_tol):
    """
    Find the matching fort.36 entry for a given CSV entry.

    Both tolerances act as hard filters:
      - freq_tol: absolute wavenumber tolerance (cm^-1) on column 10
      - int_tol: relative intensity tolerance on column 12

    Among the candidates that pass both filters, the one with the smallest
    combined relative frequency + intensity difference is selected.
    If no candidate passes both filters, no match is returned.

    Returns: (matched_cols, col10, col12) or (None, None, None)
    """
    # Convert symmetry to column value: 'e' -> 0, 'f' -> 1
    col1_sym = 0 if symmetry == 'e' else 1
    
    key = (col1_sym, J_final, J_initial)
    
    if key not in fort36_index:
        return None, None, None
    
    candidates = fort36_index[key]
    
    # Hard filter 1: frequency (absolute difference in wavenumbers)
    freq_valid = []
    for col10, col12, cols in candidates:
        if abs(col10 - frequency) <= freq_tol:
            freq_valid.append((col10, col12, cols))
    
    if not freq_valid:
        nearest = min(candidates, key=lambda c: abs(c[0] - frequency))
        print(f"  Warning: No fort.36 candidate for {symmetry}{J_initial} within "
              f"{freq_tol:.4f} cm^-1 of {frequency:.4f} cm^-1. "
              f"Nearest candidate is at {nearest[0]:.4f} cm^-1 "
              f"(offset {abs(nearest[0]-frequency):.4f}), int={nearest[1]:.3E}.")
        return None, None, None
    
    # Hard filter 2: intensity (relative difference)
    valid = []
    for col10, col12, cols in freq_valid:
        rel_int = relative_intensity_diff(col12, intensity)
        if rel_int <= int_tol:
            valid.append((col10, col12, cols, rel_int))
    
    if not valid:
        best = min(freq_valid, key=lambda c: abs(c[0] - frequency))
        print(f"  Warning: All frequency-valid candidates for {symmetry}{J_initial} "
              f"fail intensity tolerance {int_tol:.2f}. "
              f"Best candidate: col10={best[0]:.4f}, col12={best[1]:.3E} "
              f"vs experimental intensity {intensity:.3E}.")
        return None, None, None
    
    # Score among valid candidates: smallest combined relative difference wins
    best_match = None
    best_score = float('inf')
    
    for col10, col12, cols, rel_int in valid:
        rel_freq = abs(col10 - frequency) / frequency if frequency != 0 else 0.0
        score = rel_freq + rel_int
        if score < best_score:
            best_score = score
            best_match = (cols, col10, col12)
    
    return best_match


def create_comparison_table(csv_data, fort36_index, freq_tol, int_tol,
                            lev_index=None, eps=0.05, unresolved=None):
    """
    Create a comparison table matching CSV entries with fort.36 data.

    When `lev_index` is provided (i.e. the --rel_weights precision labelling is on),
    each successfully matched transition is additionally labelled by the precision
    (HP/LP) of its initial and final states. Transitions that match a fort.36 row but
    whose state pair cannot be resolved to experimental levels within `eps` are
    labelled 'Unresolved' and appended to `unresolved` (for warnings).

    Returns: list of (label, csv_freq, csv_int, fort36_freq, fort36_int, rel_diff, prec_label)
    where prec_label is None when the precision labelling is off.
    """
    results = []

    for label, csv_freq, csv_int in csv_data:
        try:
            branch, J_initial, symmetry = parse_label(label)
        except ValueError as e:
            print(f"Error parsing label '{label}': {e}")
            results.append((label, csv_freq, csv_int, None, None, None, None))
            continue

        J_final = calculate_J_final(branch, J_initial)

        matched_cols, fort36_freq, fort36_int = find_matching_entry(
            J_initial, J_final, symmetry, csv_freq, csv_int,
            fort36_index, freq_tol, int_tol
        )

        if matched_cols is None:
            print(f"No match found for {label} (J_initial={J_initial}, J_final={J_final}, sym={symmetry})")
            results.append((label, csv_freq, csv_int, None, None, None, None))
            continue

        # Calculate relative difference: (I_exp / I_calc - 1) * 100
        if fort36_int != 0:
            rel_diff = (csv_int / fort36_int - 1) * 100
        else:
            rel_diff = float('inf')

        # Optional precision labelling of the initial/final states
        prec_label = None
        if lev_index is not None:
            prec_label, details = _match_precision(matched_cols, lev_index, eps)
            if details.get('fin') is None or details.get('ini') is None:
                if unresolved is not None:
                    unresolved.append((label, csv_freq, details))

        results.append((label, csv_freq, csv_int, fort36_freq, fort36_int, rel_diff, prec_label))

    return results


def format_scientific(value, precision=3):
    """Format a number in scientific notation."""
    if value is None:
        return "N/A"
    return f"{value:.{precision}E}"


def print_comparison_table(results, show_prec=False):
    """Print the comparison table. When `show_prec` is True, a trailing
    'Precision type (Initial/Final)' column is appended."""
    prec_header = f" {'Precision type (Initial/Final)':<30}" if show_prec else ''
    width = 100 + (31 if show_prec else 0)

    # Header
    print("\n" + "=" * width)
    print("COMPARISON TABLE: Experimental vs Calculated Intensities")
    print("=" * width)
    print(f"{'Label':<10} {'Freq (exp)':<15} {'Freq (calc)':<15} {'Freq diff':<12} "
          f"{'Int (exp)':<15} {'Int (calc)':<15} {'Rel Diff (%)':<15}" + prec_header)
    print("-" * width)

    for row in results:
        label, csv_freq, csv_int = row[0], row[1], row[2]
        fort36_freq, fort36_int, rel_diff = row[3], row[4], row[5]
        prec_label = row[6] if len(row) > 6 else None
        prec_str = prec_label if prec_label is not None else 'N/A'
        prec_field = f" {prec_str:<30}" if show_prec else ''

        if fort36_freq is None or fort36_int is None:
            print(f"{label:<10} {csv_freq:<15.6f} {'N/A':<15} {'N/A':<12} "
                  f"{format_scientific(csv_int):<15} {'N/A':<15} {'N/A':<15}" + prec_field)
            continue

        freq_diff = csv_freq - fort36_freq

        print(f"{label:<10} {csv_freq:<15.6f} {fort36_freq:<15.6f} {freq_diff:<12.6f} "
              f"{format_scientific(csv_int):<15} {format_scientific(fort36_int):<15} "
              f"{rel_diff:<15.4f}" + prec_field)

    print("=" * width)


def write_table_file(results, show_prec, output_file):
    """
    Write the comparison table to `output_file` in tab-separated form (one row
    per line, columns separated by a single tab, with a header row). This layout
    can be pasted directly into a spreadsheet such as Excel, which auto-detects
    the tab as the column separator.
    """
    header = ['Label', 'Freq (exp)', 'Freq (calc)', 'Freq diff',
              'Int (exp)', 'Int (calc)', 'Rel Diff (%)']
    if show_prec:
        header.append('Precision type (Initial/Final)')

    lines = ['\t'.join(header)]
    for row in results:
        label, csv_freq, csv_int = row[0], row[1], row[2]
        fort36_freq, fort36_int, rel_diff = row[3], row[4], row[5]
        prec_label = row[6][-6:-1] if len(row) > 6 else None

        if fort36_freq is None or fort36_int is None:
            fields = [label, f'{csv_freq:.6f}', 'N/A', 'N/A',
                      format_scientific(csv_int), 'N/A', 'N/A']
        else:
            rel_str = f'{rel_diff:.4f}' if rel_diff != float('inf') else 'N/A'
            fields = [label, f'{csv_freq:.6f}', f'{fort36_freq:.6f}',
                      f'{csv_freq - fort36_freq:.6f}',
                      format_scientific(csv_int), format_scientific(fort36_int),
                      rel_str]

        if show_prec:
            fields.append(prec_label if prec_label is not None else 'N/A')
        lines.append('\t'.join(fields))

    with open(output_file, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description='Compare experimental and calculated intensities')
    parser.add_argument('--exp_file', type=str, default='input/Intensities_exp.csv',
                        help='Path to experimental CSV file (default: input/Intensities_exp.csv)')
    parser.add_argument('--calc_file', type=str, default='input/fort.36',
                        help='Path to calculated fort.36-like file (default: input/fort.36)')
    parser.add_argument('--freq_tol', type=float, default=0.01,
                        help='Absolute frequency tolerance in cm^-1 (default: 0.01)')
    parser.add_argument('--int_tol', type=float, default=0.5,
                        help='Relative intensity tolerance (default: 0.5)')
    parser.add_argument('--rel_weights', action='store_true',
                        help='label each transition by the precision (HP/LP) of its initial and '
                             'final states, using the experimental energy-level file and the '
                             'sdcalcs API. Adds a precision column, per-label counts and per-label '
                             'sd values to the output, and warns about unresolved transitions.')
    parser.add_argument('--eps', type=float, default=0.05,
                        help='tolerance (cm^-1) for matching a calculated state energy to an '
                             'experimental level; only used with --rel_weights (default: 0.05)')
    parser.add_argument('--energy_file', type=str, default='input/EnergyLevels_100_n2o.txt',
                        help='path to the experimental energy-level file with uncertainties '
                             '(N2O layout: v1 l v2 sym J E unc #lines); only used with '
                             '--rel_weights (default: input/EnergyLevels_100_n2o.txt)')
    parser.add_argument('--out_file', type=str, default=None,
                        help='if given, write the comparison table to this file as '
                             'tab-separated values (with a header row), which can be pasted '
                             'directly into a spreadsheet such as Excel')
    args = parser.parse_args()

    csv_file = args.exp_file
    fort36_file = args.calc_file

    print(f"Reading CSV file: {csv_file}")
    csv_data = parse_csv(csv_file)
    print(f"  Found {len(csv_data)} entries")

    print(f"Reading fort.36 file: {fort36_file}")
    fort36_index = parse_fort36(fort36_file)
    print(f"  Indexed {len(fort36_index)} unique quantum number combinations")

    print(f"\nFrequency tolerance: {args.freq_tol}")
    print(f"Intensity tolerance: {args.int_tol}")

    # Optional precision labelling of the initial/final states
    lev_index = None
    unresolved = []
    if args.rel_weights:
        if not os.path.isfile(args.energy_file):
            print(f"ERROR: --rel_weights is on but the energy-level file "
                  f"'{args.energy_file}' does not exist.")
            sys.exit(1)
        print(f"\nPrecision labelling: ON (state-energy tolerance eps = {args.eps} cm^-1)")
        print(f"Reading energy-level file: {args.energy_file}")
        try:
            lev_index = build_level_index(args.energy_file, dE_th=0.001)
        except ImportError:
            print("ERROR: could not import the 'sdcalcs' package. Expected it in the "
                  "sibling directory '../sd_calcs' (relative to this script).")
            sys.exit(1)
        n_levels = sum(len(v) for v in lev_index.values())
        print(f"  Indexed {n_levels} experimental levels in {len(lev_index)} (J, symmetry) groups")

    print("\nCreating comparison table...")
    results = create_comparison_table(csv_data, fort36_index, args.freq_tol, args.int_tol,
                                      lev_index=lev_index, eps=args.eps, unresolved=unresolved)

    print_comparison_table(results, show_prec=args.rel_weights)

    # Optionally write the table to a file (tab-separated, Excel-friendly)
    if args.out_file:
        write_table_file(results, args.rel_weights, args.out_file)
        print(f"\nComparison table written to: {args.out_file}")

    # Summary statistics
    valid_results = [(r[0], r[5]) for r in results if r[5] is not None and r[5] != float('inf')]
    if valid_results:
        diffs = [abs(d) for _, d in valid_results]
        avg_diff = sum(diffs) / len(diffs)
        sd = math.sqrt(sum(d ** 2 for d in diffs) / len(diffs))
        max_diff = max(diffs)
        min_diff = min(diffs)

        print(f"\nSummary:")
        print(f"  Successfully matched: {len(valid_results)}/{len(results)} entries")
        print(f"  Average absolute relative difference: {avg_diff:.4f}%")
        print(f"  Standard deviation: {sd:.4f}%")
        print(f"  Min absolute relative difference: {min_diff:.4f}%")
        print(f"  Max absolute relative difference: {max_diff:.4f}%")

        # Per-label breakdown (only when the precision labelling is on). The four
        # resolved labels are reported individually; all 'Unresolved (…)' variants
        # are grouped together under 'Unresolved'.
        if args.rel_weights:
            print(f"\nPrecision label breakdown (Initial/Final):")
            for lbl in ['HP (HP/HP)', 'LP (LP/LP)', 'Mixed (HP/LP)', 'Mixed (LP/HP)']:
                _print_label_stats(lbl, [r for r in results if r[6] == lbl])
            unres = [r for r in results if r[6] is not None and r[6].startswith('Unresolved')]
            _print_label_stats('Unresolved', unres)

    # Warnings for transitions that matched fort.36 but could not be resolved to a
    # pair of experimental levels within the state-energy tolerance.
    if args.rel_weights and unresolved:
        print(f"\nWARNING: {len(unresolved)} transition(s) could not be resolved to a pair of "
              f"experimental levels (no level within {args.eps:.3f} cm^-1 of the calculated "
              f"state energy):")
        for label, freq, det in unresolved:
            missing = []
            if det.get('fin') is None:
                missing.append(f"final state (J={det['J_final']}, sym={det['sym_final']}, "
                               f"E_calc={det['E_final']:.4f} cm^-1)")
            if det.get('ini') is None:
                missing.append(f"initial state (J={det['J_initial']}, sym={det['sym_initial']}, "
                               f"E_calc={det['E_initial']:.4f} cm^-1)")
            print(f"  {label} (freq={freq:.4f} cm^-1): no experimental "
                  + " or ".join(missing))


if __name__ == '__main__':
    main()
