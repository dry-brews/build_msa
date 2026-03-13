#!/bin/bash

### Search for sequences with homology to your query sequence via Delta-BLAST ###
### and return the filtered results as an unaligned fasta file                ###

set -e
set -o pipefail

query=./PDZ3_wt.fasta
raw_hits=./hits.txt
max_seqs=1000


# search for sequences 
./deltablast\
 -query ${query}\
 -out ${raw_hits}\
 -db nr\
 -max_target_seqs ${max_seqs}\
 -outfmt "6 qseqid sseqid qseq sseq"\
 -remote  

# filter on alignment length and write to fasta
python filter_hits.py ${raw_hits} > hits.fasta

# make sure the query sequence is included
cat ${query} hits.fasta > hits_query.fasta

# remove identical sequences
python strip_identical_seqs.py hits_query.fasta > hits_unique.fasta
