#!/usr/bin/env bash

# ------------------------------------------------------------
# This script was generated with the assistance of ChatGPT
# ------------------------------------------------------------

set -euo pipefail

# ==============================
# Usage check
# ==============================
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
THREADS=2
PARALLEL_JOBS=8

# nucmer
MINMATCH=100
CLUSTER=500
BREAKLEN=500

# delta-filter
DF_IDENTITY=90
DF_MINLEN=1000

# ==============================
# Output directories
# ==============================
OUTDIR="nucmer_by_chr"
FASTA_DIR="${OUTDIR}/fasta"
DELTA_DIR="${OUTDIR}/delta"
DF_DIR="${OUTDIR}/delta_filter"
SYRI_DIR="${OUTDIR}/syri"
LOG_DIR="${OUTDIR}/log"

mkdir -p "$FASTA_DIR" "$DELTA_DIR" "$DF_DIR" "$SYRI_DIR" "$LOG_DIR"

# ==============================
# Dependency check
# ==============================
for cmd in samtools nucmer delta-filter show-coords syri parallel; do
    if ! command -v $cmd &>/dev/null; then
        echo "ERROR: $cmd not found in PATH"
        exit 1
    fi
done

sed -i 's/\r//' "$PAIR_FILE"

[ ! -f "${REF_FA}.fai" ] && samtools faidx "$REF_FA"
[ ! -f "${QRY_FA}.fai" ] && samtools faidx "$QRY_FA"

# ==============================
# Function: run nucmer
# ==============================
run_one_pair() {
    REF_CHR=$1
    QRY_CHR=$2

    REF_SUB="${FASTA_DIR}/${REF_CHR}.fa"
    QRY_SUB="${FASTA_DIR}/${QRY_CHR}.fa"
    PREFIX="${DELTA_DIR}/${REF_CHR}_vs_${QRY_CHR}"

    samtools faidx "$REF_FA" "$REF_CHR" &>/dev/null || return
    samtools faidx "$QRY_FA" "$QRY_CHR" &>/dev/null || return

    [ ! -f "$REF_SUB" ] && samtools faidx "$REF_FA" "$REF_CHR" > "$REF_SUB"
    [ ! -f "$QRY_SUB" ] && samtools faidx "$QRY_FA" "$QRY_CHR" > "$QRY_SUB"

    nucmer \
        --maxmatch \
        -c "$CLUSTER" \
        -b "$BREAKLEN" \
        -l "$MINMATCH" \
        -t "$THREADS" \
        -p "$PREFIX" \
        "$REF_SUB" "$QRY_SUB" \
        &> "${LOG_DIR}/${REF_CHR}_vs_${QRY_CHR}.nucmer.log"
}

export -f run_one_pair
export REF_FA QRY_FA FASTA_DIR DELTA_DIR LOG_DIR
export THREADS MINMATCH CLUSTER BREAKLEN

# ==============================
# Step 1: nucmer
# ==============================
echo ">>> Running nucmer ..."
grep -v '^#' "$PAIR_FILE" | grep -v '^$' | \
parallel -j "$PARALLEL_JOBS" --colsep '\s+' \
    run_one_pair {1} {2}

cat ${DELTA_DIR}/*.delta > "${OUTDIR}/all.delta"

# ==============================
# Function: delta-filter + show-coords
# ==============================
run_delta_filter() {
    REF_CHR=$1
    QRY_CHR=$2

    BASE="${REF_CHR}_vs_${QRY_CHR}"
    DELTA_IN="${DELTA_DIR}/${BASE}.delta"
    DELTA_OUT="${DF_DIR}/${BASE}.filtered.delta"
    COORDS_OUT="${DF_DIR}/${BASE}.filtered.coords"

    [ ! -s "$DELTA_IN" ] && return

    delta-filter -1 -i "$DF_IDENTITY" -l "$DF_MINLEN" \
        "$DELTA_IN" > "$DELTA_OUT"

    show-coords -THrd "$DELTA_OUT" > "$COORDS_OUT"
}

export -f run_delta_filter
export DELTA_DIR DF_DIR DF_IDENTITY DF_MINLEN

# ==============================
# Step 2: delta-filter + show-coords
# ==============================
echo ">>> Running delta-filter + show-coords ..."
grep -v '^#' "$PAIR_FILE" | grep -v '^$' | \
parallel -j "$PARALLEL_JOBS" --colsep '\s+' \
    run_delta_filter {1} {2}

# ==============================
# Function: run syri
# ==============================
run_syri() {
    REF_CHR=$1
    QRY_CHR=$2

    BASE="${REF_CHR}_vs_${QRY_CHR}"

    COORDS="${DF_DIR}/${BASE}.filtered.coords"
    DELTA="${DF_DIR}/${BASE}.filtered.delta"
    REF_FA_SUB="${FASTA_DIR}/${REF_CHR}.fa"
    QRY_FA_SUB="${FASTA_DIR}/${QRY_CHR}.fa"

    [ ! -s "$COORDS" ] && return

    syri \
        -c "$COORDS" \
        -d "$DELTA" \
        -r "$REF_FA_SUB" \
        -q "$QRY_FA_SUB" \
        --prefix "$BASE" \
		--dir "$SYRI_DIR"
        &> "${LOG_DIR}/${BASE}.syri.log"
}

export -f run_syri
export DF_DIR SYRI_DIR FASTA_DIR

# ==============================
# Step 3: syri
# ==============================
echo ">>> Running syri ..."
grep -v '^#' "$PAIR_FILE" | grep -v '^$' | \
parallel -j "$PARALLEL_JOBS" --colsep '\s+' \
    run_syri {1} {2}

echo ">>> All done."

