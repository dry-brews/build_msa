#!/usr/bin/env bash
# jackhmmer_ebi.sh -- run jackhmmer on EBI's HMMER server and retrieve
# the aligned portions of all hit sequences as unaligned FASTA.
#
# Usage:
#   ./jackhmmer_ebi.sh <query.fasta> [output.fasta]
#
# Dependencies: curl, jq (>=1.6)

set -euo pipefail

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
API_BASE="https://www.ebi.ac.uk/Tools/hmmer/api/v1"
DATABASE="pdb"
ITERATIONS=5
POLL_INTERVAL=60    # seconds between status checks
MAX_WAIT=21600      # give up after 6 hours

# ---------------------------------------------------------------------------
# Argument handling
# ---------------------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <query.fasta> [output.fasta]" >&2
    exit 1
fi

QUERY_FILE="$1"
OUTPUT_FILE="${2:-jackhmmer_hits.fasta}"

if [[ ! -f "$QUERY_FILE" ]]; then
    echo "Error: query file '$QUERY_FILE' not found." >&2
    exit 1
fi

for cmd in curl jq; do
    command -v "$cmd" >/dev/null 2>&1 || { echo "Error: '$cmd' is required but not found." >&2; exit 1; }
done

# ---------------------------------------------------------------------------
# Read the query sequence (whole file, including header)
# ---------------------------------------------------------------------------
QUERY_SEQ=$(cat "$QUERY_FILE")

echo "==> Submitting jackhmmer job to EBI (database: $DATABASE, iterations: $ITERATIONS)..."

# ---------------------------------------------------------------------------
# Submit the job
# ---------------------------------------------------------------------------
SUBMIT_RESPONSE=$(curl -s -X POST \
    -H "Content-Type: application/json" \
    -H "Accept: application/json" \
    -d "$(jq -n \
            --arg seq "$QUERY_SEQ" \
            --arg db  "$DATABASE" \
            --argjson iters "$ITERATIONS" \
            '{input: $seq, database: $db, iterations: $iters}')" \
    "${API_BASE}/search/jackhmmer")

PARENT_JOB_ID=$(echo "$SUBMIT_RESPONSE" | jq -r '.job_id // .id // empty')

if [[ -z "$PARENT_JOB_ID" ]]; then
    echo "Error: failed to obtain a job ID. Server response:" >&2
    echo "$SUBMIT_RESPONSE" >&2
    exit 1
fi

echo "    Job ID: $PARENT_JOB_ID"
echo "    Results will be browseable at:"
echo "    https://www.ebi.ac.uk/Tools/hmmer/results/${PARENT_JOB_ID}/score"

# ---------------------------------------------------------------------------
# Poll api/v1/result/{parent_id} until all iterations are SUCCESS.
# Returns an array: [ { "id": "<iter_id>", "status": "SUCCESS", "iteration": N, ... } ]
# ---------------------------------------------------------------------------
echo "==> Waiting for iterations to complete (polling every ${POLL_INTERVAL}s)..."
ELAPSED=0
ITER_ARRAY=""

while true; do
    ITER_ARRAY=$(curl -s -H "Accept: application/json" \
        "${API_BASE}/result/${PARENT_JOB_ID}")

    IS_ARRAY=$(echo "$ITER_ARRAY" | jq 'type == "array"' 2>/dev/null || echo "false")
    if [[ "$IS_ARRAY" != "true" ]]; then
        echo "    Waiting for job to start..."
        sleep "$POLL_INTERVAL"
        ELAPSED=$(( ELAPSED + POLL_INTERVAL ))
        continue
    fi

    TOTAL_ITERS=$(echo "$ITER_ARRAY" | jq 'length')
    SUCCESS_ITERS=$(echo "$ITER_ARRAY" | jq '[.[] | select(.status == "SUCCESS")] | length')
    FAILED_ITERS=$(echo "$ITER_ARRAY" | jq '[.[] | select(.status == "FAILURE")] | length')
    PENDING_ITERS=$(echo "$ITER_ARRAY" | jq '[.[] | select(.status != "SUCCESS" and .status != "FAILURE")] | length')

    echo "    Iterations: ${SUCCESS_ITERS}/${ITERATIONS} complete..."

    if [[ "$FAILED_ITERS" -gt 0 ]]; then
        FAILED_ITER=$(echo "$ITER_ARRAY" | jq -r '[.[] | select(.status == "FAILURE")] | first | .iteration // "?"')
        echo "Error: iteration $FAILED_ITER failed on the server." >&2
        echo "    Check: https://www.ebi.ac.uk/Tools/hmmer/results/${PARENT_JOB_ID}/score" >&2
        exit 1
    fi

    if [[ "$PENDING_ITERS" -eq 0 && "$TOTAL_ITERS" -gt 0 ]]; then
        if [[ "$TOTAL_ITERS" -ge "$ITERATIONS" ]]; then
            echo "    All $ITERATIONS iterations complete."
            break
        fi
        LAST_GAINED=$(echo "$ITER_ARRAY" | jq '.[-1].convergence_stats.gained // 1')
        if [[ "$LAST_GAINED" -eq 0 ]]; then
            echo "    Converged after $TOTAL_ITERS iteration(s)."
            break
        fi
    fi

    sleep "$POLL_INTERVAL"
    ELAPSED=$(( ELAPSED + POLL_INTERVAL ))
    if [[ $ELAPSED -ge $MAX_WAIT ]]; then
        echo "Error: timed out after ${MAX_WAIT}s." >&2
        exit 1
    fi
done

# The final successful iteration's job ID holds the hit and domain results
LAST_ITER_ID=$(echo "$ITER_ARRAY" | jq -r '[.[] | select(.status == "SUCCESS")] | last | .id')
echo "    Using results from final iteration (job: $LAST_ITER_ID)"

# ---------------------------------------------------------------------------
# Fetch hit list to get accessions/descriptions, then fetch domains for
# alignment sequences.
#
# API endpoints:
#   GET api/v1/result/{iter_id}          -> { result: { hits: [...] } }
#   GET api/v1/result/{iter_id}/domains  -> { domains: [ { alignment_display: {
#                                              aseq, sqfrom, sqto, sqname, sqacc,
#                                              sqdesc, ... } } ] }
#
# The domains endpoint returns one entry per domain across all hits, in the
# same order as the hits. We use sqname/sqacc/sqfrom/sqto from alignment_display
# directly -- no need to cross-reference hits separately.
# ---------------------------------------------------------------------------
echo "==> Fetching hits with alignment data..."

# Fetch all pages of hits (paginated via ?page=N), with domain alignments embedded.
# Using /result/{iter_id}?with_domains=true returns all reported hits; the
# /domains sub-endpoint is per-hit and only returns is_included domains.
TMP_FASTA=$(mktemp /tmp/jackhmmer_hits_XXXXXX.fasta)
trap 'rm -f "$TMP_FASTA"' EXIT

TOTAL_HITS=0
PAGE=1

while true; do
    HITS_JSON=$(curl -s -H "Accept: application/json" \
        "${API_BASE}/result/${LAST_ITER_ID}?page=${PAGE}&page_size=100&with_domains=true")

    # Extract metadata in one jq call
    STATUS=$(echo "$HITS_JSON" | jq -r '.status // "unknown"')
    if [[ "$STATUS" != "SUCCESS" ]]; then
        echo "Warning: unexpected status '$STATUS' on page $PAGE" >&2
        break
    fi

    COUNT=$(echo "$HITS_JSON" | jq '.result.hits | length')
    if [[ "$COUNT" -eq 0 ]]; then
        break
    fi

    TOTAL_PAGES=$(echo "$HITS_JSON" | jq '.page_count // 1')

    # Generate all FASTA records for this page in one jq call (avoids per-hit process spawns)
    echo "$HITS_JSON" | jq -r '
        (.result.hits // [])[] |
        (if (.acc != null and .acc != "" and .acc != "null") then .acc else .name end) as $id |
        (.desc // "") as $desc |
        (.domains // [])[] |
        select(.alignment_display != null) |
        .alignment_display |
        select(.aseq != null and .aseq != "") |
        (">" + $id + "/" + (.sqfrom | tostring) + "-" + (.sqto | tostring) + " " + $desc),
        (.aseq | gsub("[.-]"; ""))
    ' >> "$TMP_FASTA"

    TOTAL_HITS=$(( TOTAL_HITS + COUNT ))
    echo "    Page $PAGE: $COUNT hits (running total: $TOTAL_HITS)"

    if [[ $PAGE -ge $TOTAL_PAGES ]]; then
        break
    fi

    PAGE=$(( PAGE + 1 ))
done

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
if [[ $TOTAL_HITS -eq 0 ]]; then
    echo "Warning: no hits found." >&2
    echo "    Check: https://www.ebi.ac.uk/Tools/hmmer/results/${PARENT_JOB_ID}/score"
    exit 0
fi

cp "$TMP_FASTA" "$OUTPUT_FILE"

echo "==> Done."
echo "    Total hit sequences : $TOTAL_HITS"
echo "    Output written to   : $OUTPUT_FILE"