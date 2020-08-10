import sys
from typing import Dict, Set, List
from glob import glob
from csv import reader as csv_reader

bowtie2_idx_extensions=["rev.1.bt2", "rev.2.bt2", "1.bt2", "2.bt2", "3.bt2", "4.bt2"]

def get_assembly(wildcards):
    assembly = glob(f'{config["assemblies_dir"]}/{wildcards.sample}*.fasta.gz')
    if len(assembly) != 1:
        raise ValueError(f'Error: expected one assembly for sample {wildcards.sample} in {config["assemblies_dir"]}, found {assembly}')
    return assembly


def get_reads(wildcards) -> List[str]:
    reads_dir=f'{config["ilmn_reads_dir"]}/{wildcards.sample}'
    reads_files = glob(f'{reads_dir}/**/*.fastq.gz')
    reads_files += glob(f'{reads_dir}/*.fastq.gz')
    if len(reads_files) == 0:
        raise FileNotFoundError(f"No reads files found in {reads_dir}")
    for read_file in reads_files:
        if " " in read_file:
            raise ValueError(f'file {read_file} has whitespace in it, this breaks the pipeline. rename the file or change the separator in the pipeline at {sys.argv[0]}')
    return reads_files


def get_samples(sample_tsv: str, excluded: Set[str], included: Set[str] = {}):
    samples=list()
    with open(sample_tsv) as f:
        reader = csv_reader(f, delimiter="\t")
        for line in reader:
            if line[0].startswith("#") or line[0].startswith("sample"):
                continue
            ID = line[0]
            if len(included) > 0: # Takes precedence over exclude set, for testing purposes
                if any(map(lambda include: include in ID, included)):
                    samples.append(ID)
            else:
                if all(map(lambda exclude: exclude not in ID, excluded)):
                    samples.append(ID)
    return samples

