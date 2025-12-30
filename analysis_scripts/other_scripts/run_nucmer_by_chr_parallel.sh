#!/usr/bin/env bash

# ------------------------------------------------------------
# This script was generated with the assistance of ChatGPT
# ------------------------------------------------------------

# Exit immediately if a command exits with a non-zero status,
# treat unset variables as an error, and propagate pipe errors
set -euo pipefail

# ==============================
# Usage check
# ==============================
# Required arguments:
#   1) Reference genome FASTA
#   2) Query genome FASTA
#   3) Chromosome pair file (two columns: ref_chr  qry_chr)
if [ $# -ne 3 ]; then
    echo "Usage: $0 <ref.fa> <qry.fa> <chr_pair.txt>"
    echo "chr_pair.txt: two columns -> ref_chr  qry_chr"
    exit 1
fi

REF_FA=$1
QRY_FA=$2
PAIR_FILE=$3

# ==============================
# Parameters
# ==============================
THREADS=2           # Threads per nucmer job
MINMATCH=100        # Minimum exact match length (-l)
CLUSTER=500         # Minimum cluster length (-c)
BREAKLEN=500        # Maximum gap between matches (-b)
PARALLEL_JOBS=8     # Number of chromosome pairs processed in parallel

# ==============================
# Output directories
# ==============================
OUTDIR="nucmer_by_chr"
FASTA_DIR="${OUTDIR}/fasta"    # Per-chromosome FASTA files
DELTA_DIR="${OUTDIR}/delta"    # nucmer delta outputs
LOG_DIR="${OUTDIR}/log"        # Log files for each comparison

mkdir -p "$FASTA_DIR" "$DELTA_DIR" "$LOG_DIR"

# ==============================
# Dependency check
# ==============================
for cmd in samtools nucmer parallel; do
    if ! command -v $cmd &>/dev/null; then
        echo "ERROR: $cmd not found in PATH"
        exit 1
    fi
done

# ==============================
# Normalize line endings in pair file (for Windows compatibility)
# ==============================
sed -i 's/\r//' "$PAIR_FILE"

# ==============================
# Build FASTA index if missing
# ==============================
[ ! -f "${REF_FA}.fai" ] && samtools faidx "$REF_FA"
[ ! -f "${QRY_FA}.fai" ] && samtools faidx "$QRY_FA"

# ==============================
# Function: run nucmer for one chromosome pair
# ==============================
run_one_pair() {
    REF_CHR=$1
    QRY_CHR=$2

    echo ">>> Processing: ${REF_CHR} vs ${QRY_CHR}"

    REF_SUB="${FASTA_DIR}/${REF_CHR}.fa"
    QRY_SUB="${FASTA_DIR}/${QRY_CHR}.fa"
    PREFIX="${DELTA_DIR}/${REF_CHR}_vs_${QRY_CHR}"

    # Check whether the chromosome exists in the FASTA
    samtools faidx "$REF_FA" "$REF_CHR" &>/dev/null || {
        echo "WARNING: ${REF_CHR} not found, skipped"
        return
    }
    samtools faidx "$QRY_FA" "$QRY_CHR" &>/dev/null || {
        echo "WARNING: ${QRY_CHR} not found, skipped"
        return
    }

    # Extract chromosome-specific FASTA if not already present
    [ ! -f "$REF_SUB" ] && samtools faidx "$REF_FA" "$REF_CHR" > "$REF_SUB"
    [ ! -f "$QRY_SUB" ] && samtools faidx "$QRY_FA" "$QRY_CHR" > "$QRY_SUB"

    # Run nucmer for this chromosome pair
    nucmer \
        --maxmatch \
        -c "$CLUSTER" \
        -b "$BREAKLEN" \
        -l "$MINMATCH" \
        -t "$THREADS" \
        -p "$PREFIX" \
        "$REF_SUB" \
        "$QRY_SUB" \
        &> "${LOG_DIR}/${REF_CHR}_vs_${QRY_CHR}.log"
}

# Export function and variables for GNU parallel
export -f run_one_pair
export REF_FA QRY_FA FASTA_DIR DELTA_DIR LOG_DIR
export THREADS MINMATCH CLUSTER BREAKLEN

# ==============================
# Run nucmer for all chromosome pairs in parallel
# ==============================
echo ">>> Running nucmer in parallel..."

grep -v '^#' "$PAIR_FILE" | grep -v '^$' | \
parallel -j "$PARALLEL_JOBS" --colsep '\s+' \
    run_one_pair {1} {2}

# ==============================
# Merge all delta files into one
# ==============================
echo ">>> Merging delta files..."
cat ${DELTA_DIR}/*.delta > "${OUTDIR}/all.delta"

echo ">>> Done."
echo "Delta : ${OUTDIR}/all.delta"
