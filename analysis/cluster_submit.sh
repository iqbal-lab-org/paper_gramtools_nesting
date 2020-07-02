#!/usr/bin/env bash

WORKFLOW=$1

function usage(){
	echo "usage: $0 workflow_name"
	exit 1
}

if [[ -z ${WORKFLOW} ]]; then usage; fi

bsub.py 5 "${WORKFLOW}" -o "run/logs/${WORKFLOW}.o" -e "run/logs/${WORKFLOW}.e" "snakemake -s analysis/workflow/${WORKFLOW}/Snakefile --use-singularity --profile lsf --verbose --latency-wait 25"

