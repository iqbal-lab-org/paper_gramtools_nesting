import click
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

sns.set_context("paper")


def usage():
    print(f"usage: {sys.argv[0]} stats_file output_dir\n")
    exit(1)


def get_NMs(condition):
    idx1, idx2 = 0, 0
    for els in condition.iteritems():
        if els[1] == "cortex_pers_ref":
            idx1 = els[0]
        if els[1] == "gramtools_genotype":
            idx2 = els[0]
    return stats_data.loc[idx1]["NM"] - stats_data.loc[idx2]["NM"]


def convert_names(df: pd.DataFrame) -> pd.DataFrame:
    gramtools_condition = [
        elem for elem in set(df["condition"]) if elem.startswith("gramtools")
    ][0]
    conversion_table = {
        "baseline_ref": "3D7",
        "cortex_baseline_ref": "cortex(3D7)",
        "samtools_baseline_ref": "samtools(3D7)",
        "cortex_pers_ref": "cortex(PR)",
        "samtools_pers_ref": "samtools(PR)",
        gramtools_condition: "gramtools",
    }
    for match, replacement in conversion_table.items():
        df = df.replace(match, replacement)
    return df


def make_condition_plot(stats_data: pd.DataFrame, metric: str, output_dir: Path):
    mean_metric = stats_data.groupby(["condition"])[metric].mean()
    condition_order = list(mean_metric.sort_values(ascending=False).index)
    for gene in set(stats_data["gene"]):
        plt.figure(figsize=(14.5, 10))
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
        plt.xlabel("")
        # ax.set_xticklabels(ax.get_xticklabels(), rotation=15)
        if metric == "NM":
            ax.set(ylabel="scaled edit distance")
        ax.figure.savefig(str(output_dir / f"py_{metric}_{gene}.pdf"))
        ax = None
    # If want to plot both box and swarmplot in facetgrid, use below, but this makes the data points and axes too small
    # plot = sns.FacetGrid(stats_data, col="gene", height=6, aspect=1)
    # plot.map(
    #    sns.boxplot, "condition", metric, order=condition_order,
    # )
    # plot.map(sns.swarmplot, "condition", metric, order=condition_order, color=".25")
    # plot.savefig(str(output_dir / f"{metric}.pdf"))


@click.command()
@click.argument("stats_file", type=click.Path(exists=True))
@click.argument("output_dir", type=click.Path())
def main(stats_file, output_dir):
    sns.set(font_scale=1.3)
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(exist_ok=True)

    stats_data = pd.read_table(str(stats_file), sep="\t")
    stats_data = convert_names(stats_data)

    # gene by condition plots for each metric
    for metric in ["NM", "delta_NM", "AS", "delta_AS"]:
        make_condition_plot(stats_data, metric, output_dir)

    ## Distribution of NM change between gramtools genotype and cortex_pers_ref
    # NM_changes = stats_data.groupby(["sample", "gene"])["condition"].agg(get_NMs)
    # NM_changes.to_csv(output_dir / "discov_genotype_NM_diffs.tsv", sep="\t")


if __name__ == "__main__":
    main()
