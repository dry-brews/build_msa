#!/usr/bin/env python3
"""
Compute per-sequence KDE weights from a distance matrix at a given bandwidth.

Weights are proportional to 1 / sum_j K(d_ij), rescaled so they sum to N_eff.
"""

import argparse
import csv
import sys
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute sequence weights from a distance matrix."
    )
    p.add_argument("--dist", required=True,
                   help="Distance matrix CSV file (square, with row/col headers)")
    p.add_argument("--bw", required=True, type=float,
                   help="KDE bandwidth (h)")
    p.add_argument("--weights", required=True,
                   help="Output TSV file: sequence_id <tab> weight")
    return p.parse_args()


def main():
    args = parse_args()
    h = args.bw

    print(f"Loading distance matrix: {args.dist}", file=sys.stderr)
    print(f"Bandwidth: h = {h}", file=sys.stderr)

    # Stream rows to keep memory at O(n) while computing O(n^2) kernel sums
    seq_ids = []
    densities = []

    with open(args.dist, newline="") as fh:
        reader = csv.reader(fh)
        next(reader)  # skip header row (column labels)
        for row in reader:
            seq_ids.append(row[0])
            d_row = np.array(row[1:], dtype=np.float64)
            kernel_row = np.exp(-(d_row ** 2) / (2 * h ** 2))
            densities.append(kernel_row.sum())

    n = len(seq_ids)
    print(f"  {n} sequences", file=sys.stderr)

    densities = np.array(densities)
    weights = 1.0 / densities
    weights /= weights.sum()
    N_eff = 1.0 / np.sum(weights ** 2)
    weights *= N_eff

    print(
        f"Weights: min={weights.min():.4e}, max={weights.max():.4e}, "
        f"effective N = {N_eff:.1f} (of {n})",
        file=sys.stderr,
    )

    with open(args.weights, "w", newline="") as fh:
        #fh.write("sequence_id\tweight\n")
        for sid, w in zip(seq_ids, weights):
            fh.write(f"{sid}\t{w:.6e}\n")

    print(f"Weights written to {args.weights}", file=sys.stderr)


if __name__ == "__main__":
    main()
