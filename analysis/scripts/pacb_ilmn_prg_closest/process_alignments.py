from pacb_ilmn_validation.process_alignments import load_gene_lengths

from pathlib import Path
from typing import Tuple

from pysam import AlignmentFile, AlignedSegment
import click


def get_gene_and_sample_name(sam_fname: Path) -> Tuple[str, str]:
    elements = sam_fname.stem.split("_")
    sample = elements[0]
    gene = "_".join(elements[1:])
    return gene, sample


class BestRead:
    def __init__(self, min_mapq: int = 0, min_qlen: int = 0):
        self.min_mapq = min_mapq
        self.best_read = None
        self.min_qlen = min_qlen

    def best_NM(self):
        if self.best_read is None:
            return "NA"
        try:
            return self.best_read.get_tag("NM")
        except KeyError:
            return "NA"

    def best_scaled_NM(self, gene_length: int):
        best_NM = self.best_NM()
        if best_NM == "NA":
            return "NA"
        return best_NM / gene_length

    def best_query(self):
        if self.best_read is None:
            return "NA"
        return self.best_read.qname

    def update(self, read: AlignedSegment):
        try:
            if read.qlen < self.min_qlen:
                return
            read_NM = read.get_tag("NM")
            if read.mapping_quality <= self.min_mapq:
                return
            best_NM = self.best_NM()
            if best_NM == "NA" or read_NM < self.best_NM():
                self.best_read = read
        except KeyError:
            return

    def output_line(self, gene_length: int):
        condition = f"closest_in_prg_mapq_{self.min_mapq}"
        best_scaled_NM = self.best_scaled_NM(gene_length)
        return f"{best_scaled_NM}\t{condition}\t{self.best_query()}"


@click.command()
@click.argument(
    "sam_fname",
    type=click.Path(exists=True),
)
@click.argument("input_bed", type=click.Path(exists=True))
@click.argument("output_file", type=str)
@click.option(
    "--min_qlen",
    type=int,
    help="minimum length of mapped query to be considered",
    default=0,
)
def main(sam_fname, input_bed, output_file, min_qlen):
    best_reads = [
        BestRead(min_mapq=0, min_qlen=min_qlen),
        BestRead(min_mapq=20, min_qlen=min_qlen),
        BestRead(min_mapq=40, min_qlen=min_qlen),
    ]
    samfile = AlignmentFile(sam_fname, "r")
    for read in samfile.fetch(until_eof=True):
        for best_read in best_reads:
            best_read.update(read)

    gene, sample = get_gene_and_sample_name(Path(sam_fname))
    gene_length = load_gene_lengths(Path(input_bed))[gene]

    with open(output_file, "w") as fout:
        for best_read in best_reads:
            output_line = best_read.output_line(gene_length)
            output_line = f"{gene}\t{sample}\t{output_line}\n"
            fout.write(output_line)


if __name__ == "__main__":
    main()
