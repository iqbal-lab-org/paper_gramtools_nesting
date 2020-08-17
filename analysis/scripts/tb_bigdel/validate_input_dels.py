"""
From a bed of deletion locations and samples they occur in (at column 4, comma-separated),
look for those deletions in those samples in PAF file (eg as output by minimap2/paftools).

Output: a tsv describing the deletions and if they are found in each sample.
"""

import click

from tb_bigdel.common import Deletion, Deletions, load_input_dels


@click.command()
@click.argument("called_PAF", type=click.Path(exists=True))
@click.argument("input_dels_bed", type=click.Path(exists=True))
@click.argument("output_file", type=str)
def main(called_vcf: click.Path, input_dels_bed: click.Path, output_file: str):
    input_dels: Deletions = load_input_dels(input_dels_bed)


if __name__ == "__main__":
    main()
