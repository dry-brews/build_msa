import sys

min_coverage = 0.85
seq_length = 112
aa_chars = set("ACDEFGHIKLMNPQRSTVWY")
gap_chars = set("-.")
filtered_seqs = 0

with open(sys.argv[1],'r') as tsv_in:
    for i, line in enumerate(tsv_in):
        row = line.strip().split('\t')
        hit_header = row[1]
        qseq = row[2]
        hseq = row[3]
        # strip gap characters out of hseq
        hseq_out = []
        for aa in hseq:
            if aa in aa_chars:
                hseq_out.append(aa)
            elif aa in gap_chars:
                continue
            else:
                sys.stderr.write("Identified non-gap, non-aa character %s, skipping\n" % aa)
                break
        if len(hseq_out) > (min_coverage * seq_length):
            sys.stdout.write(">%s\n%s\n" % (hit_header, ''.join(hseq_out)))
        else:
            filtered_seqs +=1
sys.stderr.write("Filtered %s seqs of %s total\n" % (filtered_seqs, i+1))