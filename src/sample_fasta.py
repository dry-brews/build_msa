import argparse
import random
import sys


def parse_fasta(path):
    entries = []
    header = None
    seq_lines = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    entries.append((header, "\n".join(seq_lines)))
                header = line[1:]
                seq_lines = []
            elif line:
                seq_lines.append(line)
        if header is not None:
            entries.append((header, "\n".join(seq_lines)))
    return entries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="infile", required=True)
    parser.add_argument("--num", type=int, required=True)
    parser.add_argument("--retain", nargs="*", default=[])
    args = parser.parse_args()

    entries = parse_fasta(args.infile)
    entry_map = {header: seq for header, seq in entries}

    retained = []
    for h in args.retain:
        if h in entry_map:
            retained.append((h, entry_map[h]))
        else:
            print(f"Warning: --retain header '{h}' not found in input", file=sys.stderr)

    retain_set = {h for h, _ in retained}
    remaining = [(h, s) for h, s in entries if h not in retain_set]

    sample_size = max(0, args.num - len(retained))
    if sample_size >= len(remaining):
        sampled = remaining
    else:
        sampled = random.sample(remaining, sample_size)

    for header, seq in retained + sampled:
        print(f">{header}")
        print(seq)


if __name__ == "__main__":
    main()
