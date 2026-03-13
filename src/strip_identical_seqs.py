### Minimal script that takes a fasta file that may contain identical sequences ###
### and returns a fasta file with all unique sequences. The header of the first ###
### sequence is retained if there are duplicates.                               ###

# assumes fasta format where first line starts with ">"

import sys

unique_seqs = set()
unique_count = 0
total_count  = 0

def format_seq(raw_seq):
    # raw seq is a list of lines, potentially with newlines and/or whitespace
    # return a string with 80 chars per line
    f_seq = [seq.strip() for seq in raw_seq] # strip whitespace
    f_seq = ''.join(f_seq) # make single string
    f_seq = [f_seq[i:min([i+80, len(f_seq)])] for i in range(0, len(f_seq), 80)]
    f_seq = '\n'.join(f_seq) # write printable string
    return f_seq

with open(sys.argv[1],'r') as fasta_in:
    seq = []
    for i, line in enumerate(fasta_in):
        if line.startswith(">"): # started new sequence
            total_count += 1
            if seq:
                seq = format_seq(seq) # format seq
            if i != 0 and seq not in unique_seqs: # check uniqueness
                unique_seqs.add(seq)
                sys.stdout.write(">%s\n%s\n" % (header, seq)) # print to stdout
                unique_count += 1
                                 
            header = line.strip('\n')[1:] # initialize new sequence
            seq = []
        else:
            seq.append(line.strip('\n'))
    
    # handle last sequence
    seq = format_seq(seq) # format seq
    if i != 0 and seq not in unique_seqs: # check uniqueness
        unique_seqs.add(seq)
        sys.stdout.write(">%s\n%s\n" % (header, seq)) # print to stdout
        total_count += 1

# report results
sys.stderr.write("Read %s total sequences\n" % total_count)
sys.stderr.write("Kept %s unique sequences\n" % unique_count)


