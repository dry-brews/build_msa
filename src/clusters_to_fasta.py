#!/usr/bin/env python3
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser(description="Output cluster representatives as FASTA, with optional retained sequences.")
    p.add_argument("--tsv", required=True, help="MMseqs2 cluster TSV (rep\\tmember)")
    p.add_argument("--fasta", required=True, help="FASTA file containing all sequences")
    p.add_argument("--retain", nargs="*", default=[], metavar="HEADER",
                   help="Headers to always include (replace their cluster rep if needed)")
    return p.parse_args()


def read_fasta(path):
    seqs = {}
    header = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(buf)
                header = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if header is not None:
        seqs[header] = "".join(buf)
    return seqs


def read_clusters(path):
    # rep -> list of members
    clusters = {}
    member_to_rep = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split("\t")
            rep, member = parts[0], parts[1]
            clusters.setdefault(rep, []).append(member)
            member_to_rep[member] = rep
    return clusters, member_to_rep


def main():
    args = parse_args()

    seqs = read_fasta(args.fasta)
    total_in = len(seqs)

    clusters, member_to_rep = read_clusters(args.tsv)

    # Determine which rep each retained sequence displaces
    # retained sequences appear first in output, in order given
    output_headers = []  # final ordered list of headers to print
    suppressed_reps = set()

    for retain in args.retain:
        if retain not in seqs:
            print(f"WARNING: retained header '{retain}' not found in FASTA", file=sys.stderr)
            continue
        output_headers.append(retain)
        # Find which cluster this member belongs to and suppress that rep
        rep = member_to_rep.get(retain)
        if rep and rep != retain:
            suppressed_reps.add(rep)

    # Add remaining reps (not suppressed, not already in output)
    retain_set = set(args.retain)
    for rep in clusters:
        if rep not in suppressed_reps and rep not in retain_set:
            output_headers.append(rep)

    for header in output_headers:
        if header not in seqs:
            print(f"WARNING: header '{header}' not found in FASTA, skipping", file=sys.stderr)
            continue
        print(f">{header}")
        print(seqs[header])

    print(f"Sequences read in: {total_in}", file=sys.stderr)
    print(f"Sequences printed: {len(output_headers)}", file=sys.stderr)


if __name__ == "__main__":
    main()
