from typing import Dict, List
from pathlib import Path
from collections import defaultdict
import json

import click
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap

from msps_dimorphism.site_regions import (
    Region,
    get_region,
    is_in_region,
    first_idx_in_region,
)

SiteIdx = int
Hapgs = List[int]
Hapg_Dict = Dict[SiteIdx, Hapgs]


def get_hapgs_one_site(site_json, num_samples: int) -> Hapgs:
    """
    If a sample has a null genotype, the return haplogroup (hapg)
    is -1 by convention.
    """
    result: Hapgs = [-1] * num_samples
    for sample_idx, gt in enumerate(site_json["GT"]):
        if gt[0] is not None:
            sample_hapg = site_json["HAPG"][sample_idx][0]
            result[sample_idx] = sample_hapg
    return result


def get_hapgs_all_sites(jvcf, region: Region) -> Hapg_Dict:
    num_sites = len(jvcf["Sites"])
    num_samples = len(jvcf["Samples"])
    lvl1_sites = set(jvcf["Lvl1_Sites"])
    result: Hapg_Dict = dict()

    site_is_nested = list()

    first_idx = first_idx_in_region(jvcf["Sites"], region)
    # If the first idx is not at lvl1, we can get a site following the first one which has a smaller POS and that is not in target region
    while first_idx not in lvl1_sites:
        first_idx += 1
    cur_idx = first_idx

    while is_in_region(jvcf["Sites"][cur_idx], region):
        result[cur_idx - first_idx] = get_hapgs_one_site(
            jvcf["Sites"][cur_idx], num_samples
        )
        if cur_idx in lvl1_sites:
            site_is_nested.append(False)
        else:
            site_is_nested.append(True)

        cur_idx += 1
        if cur_idx == num_sites:
            break

    result["nested"] = site_is_nested
    return result


def get_clustermap(df: pd.DataFrame, site_is_nested: List[bool]):
    """
    Use custom discrete colourmap and custom colour bar showing dicrete labels
    """
    num_colours = df.max().max() + 2  # + 2 to account for hapgs -1 and 0
    colourmap = sns.cubehelix_palette(
        reverse=True, n_colors=num_colours, hue=0.95, rot=0.5
    )

    site_is_nested = pd.Series(site_is_nested)
    col_colours = site_is_nested.map({True: "g", False: "b"})
    hmap = sns.clustermap(
        df,
        col_cluster=False,
        yticklabels=False,
        cmap=ListedColormap(colourmap),
        col_colors=col_colours,
    )
    cbar = hmap.ax_heatmap.collections[0].colorbar
    r = cbar.vmax - cbar.vmin
    cbar.set_ticks(
        [
            cbar.vmin + 0.5 * r / (num_colours) + r * i / (num_colours)
            for i in range(num_colours)
        ]
    )
    cbar.set_ticklabels(list(range(-1, num_colours - 1)))
    return hmap


@click.command()
@click.option(
    "--region",
    "-r",
    help="In the form 'SEG:start-end', as in samtools/bcftools",
    default=None,
    callback=get_region,
)
@click.option(
    "--partition_file",
    "-p",
    help="File containing sample ID partitions, tab-delimited, one partition per line",
    required=False,
    type=click.Path(exists=True),
)
@click.argument("jvcf_input", type=click.Path(exists=True))
@click.argument("output_prefix", type=str)
def main(
    region: Region,
    partition_file: click.Path,
    jvcf_input: click.Path,
    output_prefix: str,
):
    output_prefix = Path(output_prefix)

    with open(jvcf_input) as fin:
        jvcf = json.load(fin)

    hapgs_all_sites = get_hapgs_all_sites(jvcf, region)
    site_is_nested = hapgs_all_sites.pop("nested")
    df = pd.DataFrame(hapgs_all_sites)
    df.to_csv(f"{output_prefix}_hapgs.tsv", sep="\t", index=False)

    hmap = get_clustermap(df, site_is_nested)
    hmap.savefig(f"{output_prefix}_hmap.pdf")


if __name__ == "__main__":
    main()
