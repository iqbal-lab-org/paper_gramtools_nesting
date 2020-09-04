import sys
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pysam import AlignmentFile


def usage():
    print(f"usage: {sys.argv[0]} validation_stats prg_closest_stats output_dir\n")
    exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()

    output_dir = Path(sys.argv[3]).resolve()
    output_dir.mkdir(exist_ok=True)

    v_stats = pd.read_csv(sys.argv[1], sep="\t")
    c_stats = pd.read_csv(sys.argv[2], sep="\t")

    sns.set(font_scale=1.3)
    combined = pd.concat([v_stats, c_stats])

    diff_df_rows = []
    for _, group_df in combined.groupby(["gene", "sample"]):
        next_best = group_df[group_df["condition"] == "closest_in_prg"]["NM"].iloc[0]
        gtype_NM = group_df[group_df["condition"] == "gramtools_genotype"]["NM"].iloc[0]
        new_row = {
            "gene": group_df["gene"].iloc[0],
            "sample": group_df["sample"].iloc[0],
            "delta_NM": gtype_NM - next_best,
        }
        diff_df_rows.append(new_row)

    diff_df = pd.DataFrame(diff_df_rows)
    for gene in set(diff_df["gene"]):
        plt.figure(figsize=(9.5, 12))
        filtered = diff_df[diff_df["gene"] == gene]
        ax = sns.boxplot(
            data=filtered,
            y="delta_NM",
            color=sns.xkcd_rgb["windows blue"],
            whis=10000000,
        )
        ax = sns.swarmplot(data=filtered, y="delta_NM", color=".2")
        ax.figure.savefig(str(output_dir / f"{gene}_gmtools_delta.pdf"))
        ax = None
