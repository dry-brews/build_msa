#!/usr/bin/env python3
"""
Size-select homolog sequences based on a GMM fit to length distributions.

Outputs a FASTA to stdout with the wild-type sequence first, followed by
all homologs whose lengths fall within ±3σ of the GMM component that
best matches the wild-type length.
"""

import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture


def read_fasta(path):
    """Return ordered list of (header, seq) tuples."""
    records = []
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header = line[1:]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def fit_gmm(lengths, max_components=3):
    """Fit GMMs with 1..max_components and return the best by BIC."""
    X = np.array(lengths).reshape(-1, 1)
    best_gmm = None
    best_bic = np.inf
    for k in range(1, max_components + 1):
        gmm = GaussianMixture(n_components=k, random_state=42)
        gmm.fit(X)
        bic = gmm.bic(X)
        if bic < best_bic:
            best_bic = bic
            best_gmm = gmm
    return best_gmm


def parse_gmm(gmm):
    """Return (means, stds, weights) sorted by mean."""
    order = np.argsort(gmm.means_.ravel())
    means = gmm.means_.ravel()[order]
    stds = np.sqrt(gmm.covariances_.ravel()[order])
    weights = gmm.weights_[order]
    return means, stds, weights


def plot_and_save(lengths, gmm, means, stds, weights, wt_len,
                  wt_component, lower_bound, upper_bound, png_path):
    X = np.array(lengths).reshape(-1, 1)
    x = np.linspace(min(lengths) - 5, max(lengths) + 5, 1000).reshape(-1, 1)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(lengths, bins=50, color="k", density=True, alpha=0.6, label="Observed")
    ax.plot(x, np.exp(gmm.score_samples(x)), "r-", lw=2, label="GMM total")

    colors = plt.cm.tab10.colors
    for i, (m, s, w) in enumerate(zip(means, stds, weights)):
        pdf = w * (1 / (s * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x.ravel() - m) / s) ** 2)
        ls = "-" if i == wt_component else "--"
        ax.plot(x, pdf, ls=ls, color=colors[i % len(colors)],
                label=f"Component {i} (μ={m:.1f})")

    ax.axvline(lower_bound, color="blue", linestyle=":", lw=2,
               label=f"Lower bound ({lower_bound})")
    ax.axvline(upper_bound, color="blue", linestyle=":", lw=2,
               label=f"Upper bound ({upper_bound})")
    ax.axvline(wt_len, color="green", linestyle="-", lw=1.5,
               label=f"WT length ({wt_len})")

    ax.set_xlabel("Sequence length")
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(png_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--homologs", required=True,
                        help="FASTA file with homolog sequences")
    parser.add_argument("--wt", required=True,
                        help="FASTA file with a single wild-type sequence")
    parser.add_argument("--png", required=True,
                        help="Output path for the length distribution PNG")
    parser.add_argument("--n-sigma", type=float, default=3.0,
                        help="Number of sigma for bounds (default: 3)")
    args = parser.parse_args()

    # Read sequences
    wt_records = read_fasta(args.wt)
    if len(wt_records) != 1:
        sys.exit(f"ERROR: --wt file must contain exactly one sequence, found {len(wt_records)}")
    wt_header, wt_seq = wt_records[0]
    wt_len = len(wt_seq)

    homolog_records = read_fasta(args.homologs)
    lengths = [len(seq) for _, seq in homolog_records]

    print(f"WT length: {wt_len}", file=sys.stderr)
    print(f"Homologs loaded: {len(homolog_records)}", file=sys.stderr)

    # Fit GMM (BIC-selected, up to 3 components)
    gmm = fit_gmm(lengths)
    means, stds, weights = parse_gmm(gmm)
    best_k = gmm.n_components

    print(f"Best GMM components (BIC): {best_k}", file=sys.stderr)
    for i, (m, s, w) in enumerate(zip(means, stds, weights)):
        print(f"  Component {i}: mean={m:.1f}, std={s:.1f}, weight={w:.4f}",
              file=sys.stderr)

    # Identify WT component and compute bounds
    wt_component = int(np.argmin(np.abs(means - wt_len)))
    wt_mean = means[wt_component]
    wt_std = stds[wt_component]
    lower_bound = int(np.floor(wt_mean - args.n_sigma * wt_std))
    upper_bound = int(np.ceil(wt_mean + args.n_sigma * wt_std))

    print(f"WT component: mean={wt_mean:.2f}, std={wt_std:.2f}", file=sys.stderr)
    print(f"Bounds ({args.n_sigma}σ): [{lower_bound}, {upper_bound}]", file=sys.stderr)

    # Plot
    plot_and_save(lengths, gmm, means, stds, weights, wt_len,
                  wt_component, lower_bound, upper_bound, args.png)
    print(f"Plot saved to {args.png}", file=sys.stderr)

    # Filter and write output FASTA
    kept = [(h, s) for h, s in homolog_records if lower_bound <= len(s) <= upper_bound]
    print(f"Sequences retained: {len(kept)} / {len(homolog_records)} "
          f"({100 * len(kept) / len(homolog_records):.1f}%)", file=sys.stderr)

    # WT first, then filtered homologs
    sys.stdout.write(f">{wt_header}\n{wt_seq}\n")
    for header, seq in kept:
        sys.stdout.write(f">{header}\n{seq}\n")


if __name__ == "__main__":
    main()
