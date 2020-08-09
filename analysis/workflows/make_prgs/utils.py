import re
import csv
from typing import Dict
from pathlib import Path

sample_name_matcher = re.compile("([^_]+)_")


def load_vcf_names(fname: str, prg_name: str) -> Dict:
    i = 0
    vcfs = dict()
    with open(fname) as fin:
        for line in fin:
            fname = Path(config["vcf_prefix"]) / Path(
                line.strip()
            )
            if not fname.exists():
                print(f"Error: required file {fname} not found")
                exit(1)
            if prg_name.startswith("pf"):
                sample_name = sample_name_matcher.match(str(fname.name)).group(1)
            else:
                sample_name = f"sample_{i}"
                i += 1
            vcfs[sample_name] = str(fname)

    return vcfs


class Checkpoints:
    @classmethod
    def get_non_var_prg_names(cls, wildcards):
        result = []
        with checkpoints.make_beds.output["nonvar_bed"].open() as f:
            tsvreader = csv.reader(f, delimiter="\t")
            for row in tsvreader:
                result.append(row[4])
