### Perform basic QC checks on a large fasta of putative homologs ###
### 1. Remove sequences with non-aa characters                    ###
### 2. Remove duplicate (exact identical) sequences               ###
### 3. Convert lowercase characters to uppercase                  ###
### 4. Strip out gap characters ('-' or '.')                      ###
### 5. Strip commas out of header names                           ###

import sys

aa_chars       = set("ACDEFGHIKLMNPQRSTVWY")
valid_chars = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy-.')

def read_fasta(filename):
    sequences = {}
    with open(filename, 'r') as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    sequences[header] = ''.join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            sequences[header] = ''.join(seq_parts)
    return sequences


def main():
    invalid_chars = set()

    unique_seqs = set()
    comma_count     = 0
    coerce_count    = 0
    print_count     = 0
    failure_count   = 0
    nonunique_count = 0
    # read input fasta
    seqs = read_fasta(sys.argv[1])

    for i, (header, seq) in enumerate(seqs.items()):
        # Skip sequences already observed
        if seq in unique_seqs:
            nonunique_count +=1
            continue
        unique_seqs.add(seq)

        # Strip commas out of headers
        if ',' in header:
            header.replace(",", "")
            comma_count +=1
        
        # If all seq chars are uppercase aa's, good to print
        if set(seq) <= aa_chars:
            sys.stdout.write(">%s\n%s\n" % (header, seq))
            print_count +=1
            continue
        
        # Something went wrong, check for invalid characters
        if not set(seq) <= valid_chars:
            for c in set(seq) - valid_chars:
                invalid_chars.add(c)
            failure_count +=1
            continue

        # Coerce to a valid sequence
        seq = seq.upper()
        seq = seq.replace('.','')
        seq = seq.replace('-','')
        if set(seq) <= aa_chars:
            sys.stdout.write(">%s\n%s\n" % (header, seq))
            coerce_count +=1
            print_count +=1
            continue
        # if that didn't work, logic is broken
        sys.stderr.write('script failed on seq: %s\n' % seq)
        
    # report results
    sys.stderr.write("%s sequences read\n" % (i+1))
    sys.stderr.write("%s sequences written\n" % print_count)
    sys.stderr.write("%s nonunique sequences skipped\n" % nonunique_count)
    sys.stderr.write("%s commas stripped from headers\n" % comma_count)
    sys.stderr.write("%s sequences coerced to a valid format by degapping or uppercasing\n" % coerce_count)

    if len(invalid_chars) > 0:
        sys.stderr.write("%s sequences contained invalid characters and were skipped\n" % failure_count)
        sys.stderr.write("invalid chars observed: %s\n" % ' '.join(list(invalid_chars)))


if __name__ == "__main__":
    main()
