#!/usr/bin/env bash

### General file for creating a multiple sequence alignment ###
### suitable for SCA or DCA from a single protein query     ###

# Parameters
# MAX_SEQS [default: 100000] sample down to this MSA size
# MIN_SEQS [default: 1000] minimum number of sequences to continue analysis

# Steps:
# 1. Get a lot of sequences with Jackhmmer
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

# Requirements:
# hmmer (conda install -c bioconda hmmer)
# famsa (conda install -c famsa)
# mmseqs2 (sudo apt-get install mmseqs)

set -e
set -o pipefail

# 0. Parse input
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <query.fasta> [output.fasta]" >&2
    exit 1
fi

QUERY_FILE="$1"
OUTPUT_FILE="${2:-./output/alignment.fasta}"
LOG_FILE="${OUTPUT_FILE%.fasta}.log"
QUERY_HEADER=$(head -1 "$QUERY_FILE" | tr -d '\r' | sed 's/^>//' | awk '{print $1}')

if [[ ! -f "$QUERY_FILE" ]]; then
    echo "Error: query file '$QUERY_FILE' not found." >&2
    exit 1
fi

# Hard-code variables that should be user-configurable
MAX_SEQS=100000
MIN_SEQS=1000
DATABASE=/home/bryan/Workbench/alignments/databases/uniref90.fasta
CLEANUP=1

# Create output directories
mkdir -p ./tmp ./plots ./output

# Redirect all stderr to log file while still printing to terminal
exec 2> >(tee -a "$LOG_FILE" >&2)

# 1. Get sequences
echo "Searching for homologs"
jackhmmer -N 5 -o ./tmp/jh_info.txt -A ./tmp/raw_hits.sto "$QUERY_FILE" "$DATABASE"
esl-reformat fasta ./tmp/raw_hits.sto > ./tmp/raw_hits.fa

# 2. Pre-filter
echo "Applying pre-filter"
python ./src/pre-filter_sequences.py ./tmp/raw_hits.fa > ./tmp/filt_hits.fa

# 3. Filter against partials
echo "Filtering against partials"
python ./src/size_select.py --wt "$QUERY_FILE" --homologs ./tmp/filt_hits.fa --png ./plots/size_dist.png > ./tmp/sized_hits.fa

# 4. Cluster at 95% identity
echo "Clustering reads to 95% ID"
mmseqs createdb ./tmp/sized_hits.fa ./tmp/hits_db
mmseqs cluster ./tmp/hits_db ./tmp/clusters_db ./tmp/mmseqs_tmp --min-seq-id 0.95 -c 0.8 --cov-mode 0
mmseqs createtsv ./tmp/hits_db ./tmp/hits_db ./tmp/clusters_db ./tmp/clusters_db.tsv
python ./src/clusters_to_fasta.py --tsv ./tmp/clusters_db.tsv --fasta ./tmp/sized_hits.fa --retain "$QUERY_HEADER" > ./tmp/clust_reads.fa

# 5. Check > min-seqs
# Fill this in soon

# 6. Subsample reads
echo "Evaluating small subset of reads"
python ./src/sample_fasta.py --in ./tmp/clust_reads.fa --num 1000 > ./tmp/small_sample.fa

# 6a. Get distance matrix
famsa -v -gt upgma -dist_export -square_matrix ./tmp/small_sample.fa ./tmp/dist_mat.csv

# 6b. Choose bandwidth
echo "Evaluating KDE bandwidth"
BW=$(python ./src/choose_bandwidth.py --dist ./tmp/dist_mat.csv --png ./plots/seqweights_dist.png)

# 7. Take a large subset
echo "Sampling down to MAX_SEQS"
python ./src/sample_fasta.py --in ./tmp/clust_reads.fa --num $MAX_SEQS --retain "$QUERY_HEADER" > ./tmp/large_sample.fa

# 8. Align
echo "Aligning the large set, this may take a while..."
famsa -v -gt upgma ./tmp/large_sample.fa "$OUTPUT_FILE"

# 9. Weights
echo "Calculating weights, this may take a while..."
famsa -v -gt upgma -dist_export -square_matrix ./tmp/large_sample.fa ./tmp/full_dist_mat.csv
python ./src/make_weights.py --dist ./tmp/full_dist_mat.csv --bw $BW --weights ./output/weights.tsv

# Cleanup
if [[ $CLEANUP -eq 0 ]]; then
    echo "Cleaning up temporary files"
    rm -rf ./tmp
fi
