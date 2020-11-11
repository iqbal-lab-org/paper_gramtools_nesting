import sys
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context("paper")

def usage():
    print(
        f"usage: {sys.argv[0]} validation_stats prg_closest_stats output_dir gramtools_commit\n"
        "Plots the difference between gramtools-inferred-seq edit distance to truth and closest input in graph to assembly."
    )
    exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        usage()

    output_dir = Path(sys.argv[3]).resolve()
    output_dir.mkdir(exist_ok=True)

    v_stats = pd.read_csv(sys.argv[1], sep="\t")
    c_stats = pd.read_csv(sys.argv[2], sep="\t")

    gmtools_commit = sys.argv[4]

    dropped = set()
    for name in set(v_stats["condition"]):
        if name.startswith("gramtools") and gmtools_commit not in name:
            dropped.add(name)
    if len(dropped) > 0:
        print(f"Dropping entries for {dropped}")
        for name in dropped:
            v_stats = v_stats[v_stats.condition != name]

    v_stats = v_stats.replace(r"gramtools_.+", "gramtools", regex=True)

    sns.set(font_scale=1.3)
    combined = pd.concat([v_stats, c_stats])

    for condition in set(c_stats["condition"]):
        used_name = condition.replace("closest_in_prg_", "")
        diff_df_rows = []
        for _, group_df in combined.groupby(["gene", "sample"]):
            next_best = group_df[group_df["condition"] == condition]["NM"].iloc[0]
            gtype_NM = group_df[group_df["condition"] == "gramtools"]["NM"].iloc[0]
            new_row = {
                "gene": group_df["gene"].iloc[0],
                "sample": group_df["sample"].iloc[0],
                "delta_NM": gtype_NM - next_best,
            }
            diff_df_rows.append(new_row)

        diff_df = pd.DataFrame(diff_df_rows)
        for gene in set(diff_df["gene"]):
            plt.figure(figsize=(10, 12))
            filtered = diff_df[diff_df["gene"] == gene]
            ax = sns.boxplot(
                data=filtered,
                y="delta_NM",
                color=sns.xkcd_rgb["windows blue"],
                whis=10000000,
            )
            ax = sns.swarmplot(data=filtered, y=f"delta_NM", color=".2")
            ax.set(ylabel="scaled edit distance difference (gramtools vs closest input)")
            ax.figure.savefig(str(output_dir / f"{gene}_gmtools_delta_{used_name}.pdf"))
            ax = None
        diff_df.to_csv(str(output_dir / f"diff_NM_{used_name}"), sep="\t")
