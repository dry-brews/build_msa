#!/usr/bin/env python3
"""
Calculate the optimal KDE bandwidth from a FAMSA distance matrix using
leave-one-out cross-validation, then plot the resulting sequence weights.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(
        description="Find optimal KDE bandwidth via LOO-CV on a distance matrix."
    )
    p.add_argument("--dist", required=True,
                   help="Distance matrix CSV file from FAMSA (rows/cols = sequence IDs)")
    p.add_argument("--png", required=True,
                   help="Output PNG file for the sequence-weight histogram")
    return p.parse_args()


def main():
    args = parse_args()

    # --- Load distance matrix ---
    print(f"Loading distance matrix: {args.dist}", file=sys.stderr)
    dist_mat = pd.read_csv(args.dist, index_col=0).to_numpy()
    n = dist_mat.shape[0]
    print(f"  {n} sequences", file=sys.stderr)

    # --- LOO cross-validation to find optimal KDE bandwidth ---
    # Distance metric: d_ij = 1 - PID_ij  (d in [0,1], diagonal = 0)
    # Gaussian kernel: K_h(d) = exp(-d^2 / (2*h^2))
    # LOO log-likelihood: sum_i log( (1/(n-1)) * sum_{j!=i} (1/h)*K(d_ij/h) )
    # The -n*log(h) penalty prevents LL from trivially -> 0 as h -> inf.

    d_sq = dist_mat ** 2
    np.fill_diagonal(d_sq, np.nan)  # exclude self from LOO sums

    bandwidths = np.logspace(-2, 2, 80)
    loo_ll = np.zeros(len(bandwidths))

    print("Running LOO-CV over bandwidth grid...", file=sys.stderr)
    for k, h in enumerate(bandwidths):
        K = np.exp(-d_sq / (2 * h**2))
        row_sums = np.nansum(K, axis=1)
        row_sums = np.where(row_sums > 0, row_sums, np.nan)
        loo_ll[k] = np.nansum(np.log(row_sums / (n - 1))) - n * np.log(h)

    best_idx = np.nanargmax(loo_ll)
    best_h = bandwidths[best_idx]
    print(f"  Optimal bandwidth index: {best_idx} / {len(bandwidths) - 1}", file=sys.stderr)

    # --- Compute sequence weights: w_i ∝ 1 / sum_j K(d_ij) ---
    # Full kernel sum including self (d_ii = 0 -> K = 1)
    np.fill_diagonal(d_sq, 0.0)

    K_opt = np.exp(-d_sq / (2 * best_h**2))
    density = K_opt.sum(axis=1)
    weights = 1.0 / density
    weights /= weights.sum()
    N_eff = 1.0 / np.sum(weights**2)
    weights *= N_eff  # rescale so weights sum to N_eff

    print(
        f"Weights: min={weights.min():.4e}, max={weights.max():.4e}, "
        f"effective N = {N_eff:.1f} (of {n})",
        file=sys.stderr,
    )

    # --- Plot ---
    _, ax = plt.subplots(figsize=(7, 3))
    ax.hist(weights, bins=50, color="k")
    ax.set_xlabel("Weight")
    ax.set_ylabel("Count")
    ax.set_title("Sequence weights (inverse KDE density)")
    plt.tight_layout()
    plt.savefig(args.png, dpi=150)
    print(f"Histogram saved to {args.png}", file=sys.stderr)

    # --- Print result to stdout ---
    print(best_h)


if __name__ == "__main__":
    main()
