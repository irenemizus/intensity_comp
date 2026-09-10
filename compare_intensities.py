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
"""

import re
import argparse
from collections import defaultdict


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


def create_comparison_table(csv_data, fort36_index, freq_tol, int_tol):
    """
    Create a comparison table matching CSV entries with fort.36 data.
    
    Returns: list of (label, csv_freq, csv_int, fort36_freq, fort36_int, rel_diff)
    """
    results = []
    
    for label, csv_freq, csv_int in csv_data:
        try:
            branch, J_initial, symmetry = parse_label(label)
        except ValueError as e:
            print(f"Error parsing label '{label}': {e}")
            results.append((label, csv_freq, csv_int, None, None, None))
            continue
        
        J_final = calculate_J_final(branch, J_initial)
        
        matched_cols, fort36_freq, fort36_int = find_matching_entry(
            J_initial, J_final, symmetry, csv_freq, csv_int,
            fort36_index, freq_tol, int_tol
        )
        
        if matched_cols is None:
            print(f"No match found for {label} (J_initial={J_initial}, J_final={J_final}, sym={symmetry})")
            results.append((label, csv_freq, csv_int, None, None, None))
            continue
        
        # Calculate relative difference: (I_exp / I_calc - 1) * 100
        if fort36_int != 0:
            rel_diff = (csv_int / fort36_int - 1) * 100
        else:
            rel_diff = float('inf')
        
        results.append((label, csv_freq, csv_int, fort36_freq, fort36_int, rel_diff))
    
    return results


def format_scientific(value, precision=3):
    """Format a number in scientific notation."""
    if value is None:
        return "N/A"
    return f"{value:.{precision}E}"


def print_comparison_table(results):
    """Print the comparison table."""
    # Header
    print("\n" + "=" * 100)
    print("COMPARISON TABLE: Experimental vs Calculated Intensities")
    print("=" * 100)
    print(f"{'Label':<10} {'Freq (exp)':<15} {'Freq (calc)':<15} {'Freq diff':<12} "
          f"{'Int (exp)':<15} {'Int (calc)':<15} {'Rel Diff (%)':<15}")
    print("-" * 100)
    
    for label, csv_freq, csv_int, fort36_freq, fort36_int, rel_diff in results:
        if fort36_freq is None or fort36_int is None:
            print(f"{label:<10} {csv_freq:<15.6f} {'N/A':<15} {'N/A':<12} "
                  f"{format_scientific(csv_int):<15} {'N/A':<15} {'N/A':<15}")
            continue
        
        freq_diff = csv_freq - fort36_freq
        
        print(f"{label:<10} {csv_freq:<15.6f} {fort36_freq:<15.6f} {freq_diff:<12.6f} "
              f"{format_scientific(csv_int):<15} {format_scientific(fort36_int):<15} "
              f"{rel_diff:<15.4f}")
    
    print("=" * 100)


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
    
    print("\nCreating comparison table...")
    results = create_comparison_table(csv_data, fort36_index, args.freq_tol, args.int_tol)
    
    print_comparison_table(results)
    
    # Summary statistics
    valid_results = [(r[0], r[5]) for r in results if r[5] is not None and r[5] != float('inf')]
    if valid_results:
        diffs = [abs(d) for _, d in valid_results]
        avg_diff = sum(diffs) / len(diffs)
        max_diff = max(diffs)
        min_diff = min(diffs)
        
        print(f"\nSummary:")
        print(f"  Successfully matched: {len(valid_results)}/{len(results)} entries")
        print(f"  Average absolute relative difference: {avg_diff:.4f}%")
        print(f"  Min absolute relative difference: {min_diff:.4f}%")
        print(f"  Max absolute relative difference: {max_diff:.4f}%")


if __name__ == '__main__':
    main()
