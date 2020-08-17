"""
From a bed of deletion locations and samples they occur in (at column 4, comma-separated),
look for those deletions in those samples in a multi-sample vcf.

Output: a tsv describing the deletions and if they are found in each sample.
"""
from pathlib import Path

import click
from pysam import VariantFile

from tb_bigdel.common import Deletion, Deletions, load_input_dels


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
            if gtyped_del.overlaps(input_del):
                found_dels.append(gtyped_del)
        if len(found_dels) > 1:
            print(
                f"WARNING: input del {input_del} found in >1 separate records: {found_dels}"
            )
        found_del = found_dels[0] if len(found_dels) > 0 else None
        found_samples = set()
        for found in found_dels:
            found_samples.update(found.samples)
        if found_del is not None:
            delta_len = input_del.del_len - found_del.del_len
            delta_pos = input_del.start - found_del.start
        for sample in input_del.samples:
            line = f"{input_del.start}\t{input_del.del_len}\t{sample}\t"
            if found_del is None:
                line += "0\t.\t.\n"
            elif sample not in found_samples:
                line += f"0\t{delta_len}\t{delta_pos}\n"
            else:
                line += f"1\t{delta_len}\t{delta_pos}\n"
            fout.write(line)

    fout.close()


if __name__ == "__main__":
    main()
