#!/usr/bin/env bash
bsub.py 5 "msps_dimorphism" -o run/logs/msps_dimorphism.o -e run/logs/msps_dimorphism.e "snakemake -s analysis/workflow/msps_dimorphism/Snakefile --use-singularity --profile lsf --verbose --latency-wait 10"

