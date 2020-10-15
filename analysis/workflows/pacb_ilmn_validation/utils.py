import sys
from typing import Dict, Set, List
from glob import glob
from csv import reader as csv_reader

bowtie2_idx_extensions = ["rev.1.bt2", "rev.2.bt2", "1.bt2", "2.bt2", "3.bt2", "4.bt2"]


def get_assembly(wildcards):
    assembly = glob(f'{config["assemblies_dir"]}/{wildcards.sample}*.fasta.gz')
    if len(assembly) != 1:
        raise ValueError(
            f'Error: expected one assembly for sample {wildcards.sample} in {config["assemblies_dir"]}, found {assembly}'
        )
    return assembly


def get_samples(
    sample_tsv: str, excluded: Set[str] = set(), included: Set[str] = set()
):
    samples = list()
    with open(sample_tsv) as f:
        reader = csv_reader(f, delimiter="\t")
        for line in reader:
            if line[0].startswith("#") or line[0].startswith("sample"):
                continue
            ID = line[0]
            if (
                len(included) > 0
            ):  # Takes precedence over exclude set, for testing purposes
                if any(map(lambda include: include in ID, included)):
                    samples.append(ID)
            else:
                if all(map(lambda exclude: exclude not in ID, excluded)):
                    samples.append(ID)
    return samples
