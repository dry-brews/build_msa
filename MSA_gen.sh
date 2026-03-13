#!/usr/bin/env bash

### General file for creating a multiple sequence alignment ###
### suitable for SCA or DCA from a single protein query     ###

# Parameters
# MAX_SEQS [default: 100000] sample down to this MSA size
# MIN_SEQS [default: 1000] minimum number of sequences to continue analysis

# Steps:
# 1. Get a lot of sequences with Jackhmmer (call_jackhmmer.sh)
# 2. Pre-filter to remove duplicates and ambiguous hits (pre-filter_sequences.py)
# 3. Filter against partials, i.e., require long aligned length (size_select.py) (also adds query to set)
# 4. Cluster at 95% identity to collapse very similar sequencesm (MMseqs2)
#     a. Covert to fasta and ensure wild type is included at top
# 5. Ensure at least ${min-seqs} sequences or abort
# 6. Take a random subsample of 1k sequences (sample_fasta.py)
#     a. Align and return a complete distance matrix (FAMSA)
#     b. Choose an optimal KDE bandwidth from the subsample and return as a variable (choose_bandwidth.py)
# 7. Take a random subsample of {max-seqs} sequences (sample_fasta.py)
# 8. Align the large set (FAMSA)
# 9. Compute the large distance matrix (FAMSA) and convert to sequence weights (make_weights.py)

# 0. Parse input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <query.fasta> [output.fasta]" >&2
    exit 1
fi

QUERY_FILE="$1"
OUTPUT_FILE="${2:-alignment.fasta}"
LOG_FILE="${OUTPUT_FILE%.fasta}.log"

if [[ ! -f "$QUERY_FILE" ]]; then
    echo "Error: query file '$QUERY_FILE' not found." >&2
    exit 1
fi

# Redirect all stderr to log file while still printing to terminal
exec 2> >(tee -a "$LOG_FILE" >&2)

## Extract the sequence ID from the first header line of the query FASTA
QUERY_HEADER=$(head -1 "$QUERY_FILE" | sed 's/^>//' | awk '{print $1}')
MAX_SEQS=100000
MIN_SEQS=1000

# 1. Get sequences
echo "Searching for homologs"
./src/call_jackhmmer.sh "$QUERY_FILE" ./raw_hits.fa

# 2. Pre-filter
echo "Applying pre-filter"
python ./src/pre-filter_sequences.py ./raw_hits.fa > ./filt_hits.fa

# 3. Filter against partials
echo "Filtering against partials"
python ./src/size_select.py --wt "$QUERY_FILE" --homologs ./filt_hits.fa --png size_dist.png > sized_hits.fa

# 4. Cluster at 95% identity
echo "Clustering reads to 95% ID"
mmseqs createdb sized_hits.fa hits_db
mmseqs cluster hits_db clusters_db ./tmp --min-seq-id 0.95 -c 0.8 --cov-mode 0
mmseqs createtsv hits_db hits_db clusters_db clusters_db.tsv
python ./src/clusters_to_fasta.py --tsv clusters_db.tsv --fasta sized_hits.fa --retain "$QUERY_HEADER" > clust_reads.fa

# 5. Check > min-seqs
# Fill this in soon

# 6. Subsample reads
echo "Evaluating small subset of reads"
python ./src/sample_fasta.py --in clust_reads.fa --num 1000 > small_sample.fa

# 6a. Get distance matrix
famsa -v -gt upgma -dist_export -square_matrix small_sample.fa dist_mat.csv

# 6b. Choose bandwidth
echo "Evaluating KDE bandwidth"
BW=$(python ./src/choose_bandwidth.py --dist dist_mat.csv --png seqweights_dist.png)

# 7. Take a large subset
echo "Sampling down to MAX_SEQS"
python ./src/sample_fasta.py --in clust_reads.fa --num $MAX_SEQS --retain "$QUERY_HEADER" > large_sample.fa

# 8. Align
echo "Aligning the large set, this may take a while..."
famsa -v -gt upgma large_sample.fa $OUTPUT_FILE

# 9. Weights
echo "Calculating weights, this may take a while..."
famsa -v -gt upgma -dist_export -square_matrix large_sample.fa full_dist_mat.csv
python ./src/make_weights.py --dist full_dist_mat.csv --bw $BW --weights weights.tsv