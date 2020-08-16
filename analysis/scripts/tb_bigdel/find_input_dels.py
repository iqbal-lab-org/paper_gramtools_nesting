"""
From a bed of deletion locations and samples they occur in (at column 4, comma-separated),
look for those deletions in those samples in a multi-sample vcf.

Output: a tsv describing the deletions and if they are found in each sample.
"""
from typing import List
from pathlib import Path

import click
from pysam import VariantFile


class Deletion:
    def __init__(self, start: int, stop: int, del_len: int, samples: set):
        self.start = start
        self.stop = stop
        self.del_len = del_len
        self.samples = samples

    def modify_by(self, left_delta: int, right_delta: int):
        self.start += left_delta
        self.stop += right_delta

    def _contains(self, position: int):
        return self.start <= position <= self.stop

    def spans(self, other: "Deletion"):
        return self._contains(other.start) and self._contains(other.stop)

    def compatible_with(self, other: "Deletion"):
        return self._contains(other.start) or self._contains(other.stop)

    def __len__(self) -> int:
        return self.stop - self.start + 1

    def __lt__(self, other: "Deletion") -> bool:
        return self.start < other.start

    def __eq__(self, other: "Deletion") -> bool:
        return self.start == other.start and self.stop == other.stop

    def __repr__(self):
        return f"[{self.start}, {self.stop}]: {self.samples}"


Deletions = List[Deletion]


def load_input_dels(input_dels_bed) -> Deletions:
    input_dels: Deletions = list()
    with Path(input_dels_bed).open() as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            rows = line.split("\t")
            samples = set(rows[3].strip().split(","))
            start, stop = int(rows[1]) + 1, int(rows[2])  # +1 : bed start is 0-based
            del_len = stop - start  # That's how I picked them: len(alt) == 1
            input_dels.append(Deletion(start, stop, del_len, samples))

    return input_dels


def load_gtyped_dels(called_vcf) -> Deletions:
    gtyped_dels: Deletions = list()
    vcf_recs = VariantFile(called_vcf)
    for rec in vcf_recs.fetch():
        samples = set()
        bigdels = []
        for name, vals in rec.samples.items():
            gt_allele = vals["GT"][0]
            if gt_allele is None:
                continue
            if len(rec.ref) - len(rec.alleles[gt_allele]) > 100:
                samples.add(name)
                bigdels.append(rec.alleles[gt_allele])

        if len(samples) == 0:
            continue

        del_len = len(rec.ref) - min(map(len, bigdels))  # Take len of largest deletion

        gtyped_dels.append(
            Deletion(rec.pos, rec.pos + len(rec.ref) - 1, del_len, samples)
        )
    return gtyped_dels


@click.command()
@click.argument("called_vcf", type=click.Path(exists=True))
@click.argument("input_dels_bed", type=click.Path(exists=True))
@click.argument("output_file", type=str)
def main(called_vcf: click.Path, input_dels_bed: click.Path, output_file: str):

    input_dels: Deletions = load_input_dels(input_dels_bed)
    gtyped_dels: Deletions = load_gtyped_dels(called_vcf)

    fout = open(output_file, "w")
    for input_del in input_dels:
        found_dels = []
        for gtyped_del in gtyped_dels:
            if gtyped_del.compatible_with(input_del):
                found_dels.append(gtyped_del)
        if len(found_dels) > 1:
            raise ValueError(
                f"input del {input_del} found in >1 separate records: {found_dels}"
            )
        found_del = found_dels[0] if len(found_dels) == 1 else None
        if found_del is not None:
            delta_len = input_del.del_len - found_del.del_len
            delta_pos = input_del.start - found_del.start
        for sample in input_del.samples:
            line = f"{input_del.start}\t{input_del.del_len}\t{sample}\t"
            if found_del is None:
                line += "0\t.\t.\n"
            elif sample not in found_del.samples:
                line += f"0\t{delta_len}\t{delta_pos}\n"
            else:
                line += f"1\t{delta_len}\t{delta_pos}\n"
            fout.write(line)

    fout.close()


if __name__ == "__main__":
    main()
