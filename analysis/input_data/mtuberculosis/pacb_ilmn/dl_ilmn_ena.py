#!/usr/bin/env python3

import csv
from pathlib import Path
import shutil
import subprocess

def get_one_run(run_id, outdir):
    command = f"enaDataGet --format fastq {run_id}"
    print("Start:", command)
    subprocess.check_output(command, shell=True, cwd=outdir)
    print("Finish:", command)



def get_one_sample(sample_name, run_ids):
    outdir = Path("ilmn_reads") / sample_name
    outdir.mkdir(parents=True,exist_ok=True)

    exist = [Path(f"{outdir}/reads_{i}.fastq.gz").exists() for i in [1,2]]
    if False not in exist:
        print("{outdir} already has reads, skipped")
        return

    for run in run_ids:
        get_one_run(run, outdir)

    for i in [1,2]:
        if len(run_ids) > 1:
            reads = " ".join([f"{x}/{x}_{i}.fastq.gz" for x in run_ids])
            command = f"cat {reads} > reads_{i}.fastq.gz"
        else:
            fname=Path(f"{run_ids[0]}/{run_ids[0]}_{i}.fastq.gz")
            if fname.exists():
                command = f"mv {fname} reads_{i}.fastq.gz"
            else:
                command = f"cp {run_ids[0]}/{run_ids[0]}.fastq.gz reads.fastq.gz"

        print("Start:", command)
        subprocess.check_output(command, shell=True, cwd=outdir)
        print("Finish:", command)

    for run_id in run_ids:
        shutil.rmtree(outdir / run_id)



with open("illumina_run_ids.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for d in reader:
        runs = d["runs"].split(",")
        get_one_sample(d["sample"], runs)

