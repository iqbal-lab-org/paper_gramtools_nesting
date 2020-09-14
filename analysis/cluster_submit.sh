#!/usr/bin/env bash

set -eu

WORKFLOW=$1

function usage(){
    echo "usage: $0 workflow_name"
    exit 1
}

if [[ -z ${WORKFLOW} ]]; then usage; fi

SINGULARITY_BINDS="/hps/nobackup2/iqbal,/nfs/leia/research/iqbal"
SINGULARITY_ARGS="--contain --bind $SINGULARITY_BINDS"

LOG_DIR="analysis/logs"
MEMORY=5000

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -M "$MEMORY" \
    -o "${LOG_DIR}/${WORKFLOW}.o" \
    -e "${LOG_DIR}/${WORKFLOW}.e" \
    -J "${WORKFLOW}_snakemake" \
    snakemake -s analysis/workflows/${WORKFLOW}/Snakefile \
    --profile lsf --verbose --latency-wait 25 \
    --use-singularity --singularity-args "$SINGULARITY_ARGS"

