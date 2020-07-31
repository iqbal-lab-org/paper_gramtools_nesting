import sys
from collections import namedtuple
from typing import Dict, Tuple, List
from pathlib import Path

from pysam import AlignmentFile
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def usage():
    print(
        f"usage: {sys.argv[0]} input_dir output_dir input_bed\n"
        "The input_dir should contain the .sam files to analyse"
    )
    exit(1)


class MultipleAlignmentsError(Exception):
    pass


Scores = namedtuple("Scores", ["NM", "AS"])
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
        NM = read.get_tag("NM")
        if len(gene_lengths) > 0:
            NM = NM / gene_lengths[gene_name]
        AS = (
            read.get_tag("AS") * -1
        )  # Convert to positive so that high is bad, low is good, like for NM.
        result[gene_name] = Scores(NM, AS)
    return result


def compute_delta(scores_1: Scores, scores_2: Scores) -> Scores:
    newattrs = dict()
    for attr in Scores._fields:
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
        fieldnames = ["sample", "gene", "condition", "NM", "AS", "delta_NM", "delta_AS"]
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
                    delta_scores[gene].NM,
                    delta_scores[gene].AS,
                ]
                stats_file.write("\t".join(map(str, row)) + "\n")


def make_condition_plot(stats_data: pd.DataFrame, metric: str, output_dir: Path):
    mean_metric = stats_data.groupby(["condition"])[metric].mean()
    condition_order = list(mean_metric.sort_values(ascending=False).index)
    for gene in set(stats_data["gene"]):
        plt.figure(figsize=(10, 7))
        filtered = stats_data[stats_data["gene"] == gene]
        ax = sns.boxplot(
            data=filtered,
            x="condition",
            y=metric,
            order=condition_order,
            color=sns.xkcd_rgb["windows blue"],
            whis=10000000,
        )
        ax = sns.swarmplot(
            data=filtered, x="condition", y=metric, color=".2", order=condition_order
        )
        ax.figure.savefig(str(output_dir / f"{metric}_{gene}.pdf"))
        ax = None
    # If want to plot both box and swarmplot in facetgrid, use below, but this makes the data points and axes too small
    # plot = sns.FacetGrid(stats_data, col="gene", height=6, aspect=1)
    # plot.map(
    #    sns.boxplot, "condition", metric, order=condition_order,
    # )
    # plot.map(sns.swarmplot, "condition", metric, order=condition_order, color=".25")
    # plot.savefig(str(output_dir / f"{metric}.pdf"))


if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()

    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(exist_ok=True)
    output_stats = output_dir / "stats.tsv"

    if not output_stats.exists():
        input_dir = Path(sys.argv[1]).resolve()
        if not input_dir.exists():
            print(f"Error: {input_dir} not found")
            usage()

        sam_file_list = list(input_dir.glob(f"*.sam"))
        if len(sam_file_list) == 0:
            print(f"Error: no .sam files in {input_dir}")
            usage()

        input_bed = Path(sys.argv[3]).resolve()
        if not input_bed.exists():
            print(f"Error: {input_bed} not found")
            usage()
        gene_lengths = load_gene_lengths(input_bed)

        write_stats(sam_file_list, output_stats, gene_lengths)
    else:
        print(f"Found existing {output_stats} so re-using it- delete it to regenerate")

    ### Write plots ###
    stats_data = pd.read_table(str(output_stats), sep="\t")

    for cond in ["NM", "delta_NM", "AS", "delta_AS"]:
        make_condition_plot(stats_data, cond, output_dir)

    ## Analyse stats ##
    ## Get the distribution of NM change between gramtools genotype and cortex_pers_ref
    data = pd.read_csv(str(output_stats), sep="\t")

    def get_NMs(condition):
        idx1, idx2 = 0, 0
        for els in condition.iteritems():
            if els[1] == "cortex_pers_ref":
                idx1 = els[0]
            if els[1] == "gramtools_genotype":
                idx2 = els[0]
        return data.loc[idx1]["NM"] - data.loc[idx2]["NM"]

    NM_changes = data.groupby(["sample", "gene"])["condition"].agg(get_NMs)
    NM_changes.to_csv(output_dir / "discov_genotype_NM_diffs.tsv", sep="\t")
