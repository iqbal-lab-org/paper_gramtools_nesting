#!/usr/bin/env bash

WORKFLOW=$1

function usage(){
	echo "usage: $0 workflow_name"
	exit 1
}

if [[ -z ${WORKFLOW} ]]; then usage; fi

bsub.py 5 "${WORKFLOW}" -o "analysis/logs/${WORKFLOW}.o" -e "analysis/logs/${WORKFLOW}.e" "snakemake -s analysis/workflows/${WORKFLOW}/Snakefile --use-singularity --profile lsf --verbose --latency-wait 25"

