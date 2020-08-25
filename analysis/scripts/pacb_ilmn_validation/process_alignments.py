import sys
from collections import namedtuple
from typing import Dict, Tuple, List
from pathlib import Path

from pysam import AlignmentFile


def usage():
    print(
        f"usage: {sys.argv[0]} input_dir input_bed output_dir\n"
        "The input_dir should contain the .sam files to analyse.\n"
        "The input_bed should contain the read names in column 4."
    )
    exit(1)


class MultipleAlignmentsError(Exception):
    pass


Scores = namedtuple("Scores", ["NM", "AS", "MAPQ"])
GeneScores = Dict[str, Scores]
BaselineScores: Dict[str, GeneScores] = dict()


def get_sample_and_condition_name(sam_fname: Path) -> Tuple[str, str]:
    elements = sam_fname.stem.split("_")
    sample = elements[-1]
    condition = "_".join(elements[:-1])
    return sample, condition


def load_gene_lengths(bed_fname: Path) -> Dict[str, int]:
    result = dict()
    with bed_fname.open() as fin:
        for line in fin:
            cols = line.split("\t")
            result[cols[3].strip()] = int(cols[2]) - int(cols[1]) + 1
    return result


def get_scores(sam_fname: Path, gene_lengths: Dict[str, int]) -> GeneScores:
    result: GeneScores = dict()
    samfile = AlignmentFile(sam_fname, "r")
    for read in samfile.fetch(until_eof=True):
        gene_name = read.query_name
        if gene_name in result:
            raise MultipleAlignmentsError(
                f"Several alignments to gene {gene_name} in file {sam_fname}"
            )
        try:
            NM = read.get_tag("NM")
            NM = NM / gene_lengths[gene_name]
            # Convert to positive so that high is bad, low is good, like for NM.
            AS = read.get_tag("AS") * -1
            MAPQ = read.mapping_quality
            result[gene_name] = Scores(NM, AS, MAPQ)
        except KeyError:  # Case: read not aligned
            result[gene_name] = Scores("NA", "NA", "NA")

    return result


def compute_delta(scores_1: Scores, scores_2: Scores) -> Scores:
    newattrs = dict()
    for attr in Scores._fields:
        score_1 = getattr(scores_1, attr)
        score_2 = getattr(scores_2, attr)
        if score_1 == "NA" or score_2 == "NA":
            newattrs[attr] = "NA"
        else:
            newattrs[attr] = getattr(scores_2, attr) - getattr(scores_1, attr)
    return Scores(**newattrs)


def get_delta_scores(baseline_genes: GeneScores, query_genes: GeneScores) -> GeneScores:
    """Computes the scores of genes in `query_genes` minus the scores of genes in `baseline_genes`"""
    if baseline_genes.keys() != query_genes.keys():
        raise ValueError(
            f"Cannot get score deltas for dicts without same keys: {baseline_genes} vs {query_genes}"
        )
    result: GeneScores = dict()
    for gene_name, scores_1 in baseline_genes.items():
        scores_2 = query_genes[gene_name]
        result[gene_name] = compute_delta(scores_1, scores_2)
    return result


def write_stats(sam_file_list: List[Path], output_stats: Path, gene_lengths):

    # First pass: get baseline scores
    baseline_scores: BaselineScores = dict()
    for sam_fname in sam_file_list:
        sample, condition = get_sample_and_condition_name(sam_fname)
        if condition == "baseline_ref":
            baseline_scores[sample] = get_scores(sam_fname, gene_lengths)

    # Second pass: get scores and deltas relative to baseline
    with output_stats.open("w") as stats_file:
        fieldnames = [
            "sample",
            "gene",
            "condition",
            "NM",
            "AS",
            "MAPQ",
            "delta_NM",
            "delta_AS",
            "delta_MAPQ",
        ]
        stats_file.write("\t".join(fieldnames) + "\n")

        for sam_fname in sam_file_list:
            sample, condition = get_sample_and_condition_name(sam_fname)
            scores = get_scores(sam_fname, gene_lengths)
            delta_scores = get_delta_scores(baseline_scores[sample], scores)
            for gene in delta_scores:
                row = [
                    sample,
                    gene,
                    condition,
                    scores[gene].NM,
                    scores[gene].AS,
                    scores[gene].MAPQ,
                    delta_scores[gene].NM,
                    delta_scores[gene].AS,
                    delta_scores[gene].MAPQ,
                ]
                stats_file.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()

    output_dir = Path(sys.argv[3]).resolve()
    output_dir.mkdir(exist_ok=True)
    output_stats = output_dir / "stats.tsv"

    if not output_stats.exists():
        input_dir = Path(sys.argv[1]).resolve()
        input_bed = Path(sys.argv[2]).resolve()
        for _input in [input_dir, input_bed]:
            if not _input.exists():
                print(f"Error: {input_bed} not found")
                usage()

        sam_file_list = list(input_dir.glob(f"*.sam"))
        if len(sam_file_list) == 0:
            print(f"Error: no .sam files in {input_dir}")
            usage()

        gene_lengths = load_gene_lengths(input_bed)

        write_stats(sam_file_list, output_stats, gene_lengths)
    else:
        print(f"Found existing {output_stats}, nothing to do. Delete it to regenerate")
