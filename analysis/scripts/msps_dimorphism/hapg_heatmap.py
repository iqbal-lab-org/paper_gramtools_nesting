from typing import Dict, List
from pathlib import Path
from collections import defaultdict
import json

import click
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
from matplotlib.pyplot import gcf

from jvcf_processing import (
    Region,
    click_get_region,
    is_in_region,
    first_idx_in_region,
)

SiteIdx = int
Hapgs = List[int]
Hapg_Dict = Dict[SiteIdx, Hapgs]

country_to_colour = {"Ghana": "gold", "Laos": "firebrick", "Cambodia": "royalblue"}
nested_colour_mapping = {True: "gray", False: "navajowhite"}


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


def get_country_colouring(metadata_fname, df):
    """Loads sample metadata and assigns country colouring"""
    metadata = pd.read_csv(metadata_fname,sep="\t")
    sample_to_country = dict(zip(list(metadata["sample"]),list(metadata["country"])))
    country_colours = list()
    for sname in list(df.index):
        country_colours.append(country_to_colour[sample_to_country[sname]])
    return country_colours

def get_cell_colouring(df):
    """Custom cell colouring based on viridis palette"""
    num_colours = df.max().max() + 1  # + 1 to account for haplogroup of 0
    min_colour, max_colour = 30, 100
    to_map=list(map(lambda val: val/100,range(min_colour, max_colour, int((max_colour - min_colour)/num_colours))))
    if len(to_map) > num_colours:
        to_map = to_map[:-1]
    cmap=sns.color_palette("viridis", as_cmap=True)
    # There are actually num_colours + 1 haplogroup values as haplogroup can be -1; we colour -1 differently below
    null_hapg_colour = (0.18195582, 0.11955283, 0.23136943) # blackish colour (RGB)
    colourmap = [null_hapg_colour] + [cmap(x) for x in to_map]
    return colourmap

def add_clustermap_legends(hmap, num_colours):
    """Adds row, colour, and cell colour legends to clustermap
    See https://stackoverflow.com/q/27988846/12519542 for some of the code for this
    """
    cbar = hmap.ax_heatmap.collections[0].colorbar
    r = cbar.vmax - cbar.vmin
    cbar.set_ticks(
        [
            cbar.vmin + 0.5 * r / (num_colours) + r * i / (num_colours)
            for i in range(num_colours)
        ]
    )
    cbar.set_ticklabels(list(range(-1, num_colours - 1)))


    for label,colour in nested_colour_mapping.items():
        hmap.ax_col_dendrogram.bar(0, 0, color=colour, label=label, linewidth=0)

    l1 = hmap.ax_col_dendrogram.legend(title='Nested site', loc="center", ncol=1, bbox_to_anchor=(0.6, 0.88), bbox_transform=gcf().transFigure)

    for label,colour in country_to_colour.items():
        hmap.ax_row_dendrogram.bar(0, 0, color=colour, label=label, linewidth=0)

    l2 = hmap.ax_row_dendrogram.legend(title='Sample country', loc="center", ncol=1, bbox_to_anchor=(0.25, 0.88), bbox_transform=gcf().transFigure)


@click.command()
@click.option(
    "--region",
    "-r",
    help="In the form 'SEG:start-end', as in samtools/bcftools",
    default=None,
    callback=click_get_region,
)
@click.option(
    "--partition_file",
    "-p",
    help="File containing sample ID partitions, tab-delimited, one partition per line",
    required=False,
    type=click.Path(exists=True),
)
@click.argument("jvcf_input", type=click.Path(exists=True))
@click.argument("metadata_file", type=click.Path(exists=True))
@click.argument("output_prefix", type=str)
def main(
    region: Region,
    partition_file: click.Path,
    jvcf_input: click.Path,
    metadata_file: click.Path,
    output_prefix: str,
):
    output_prefix = Path(output_prefix)

    with open(jvcf_input) as fin:
        jvcf = json.load(fin)


    hapgs_all_sites = get_hapgs_all_sites(jvcf, region)
    site_is_nested = hapgs_all_sites.pop("nested")
    sample_names = [sample["Name"] for sample in jvcf["Samples"]]
    df = pd.DataFrame(hapgs_all_sites, index=sample_names)
    df.to_csv(f"{output_prefix}_hapgs.tsv", sep="\t")

    #### Get colourings ####
    ## Cell colouring: according to haplogroup
    cell_colourmap = get_cell_colouring(df)

    ## Row colouring: colour samples by country of origin
    country_colours = get_country_colouring(metadata_file, df)

    ## Column colouring: colour sites by whether they are nested
    site_is_nested = pd.Series(site_is_nested)
    site_colours = site_is_nested.map(nested_colour_mapping)

    ## row-clustered clustermap
    hmap_clustered = sns.clustermap(
    df,
    method="average",
    metric="euclidean",
    col_cluster=False,
    yticklabels=False,
    cmap=cell_colourmap,
    col_colors=site_colours,
    row_colors=country_colours,
    cbar_kws={"label":"haplogroup"},
    cbar_pos=(1.05, 0.35, 0.05, 0.18)
)
    add_clustermap_legends(hmap_clustered, len(cell_colourmap))
    hmap_clustered.savefig(f"{output_prefix}_hmap_clustered.pdf")

    ## non-row-clustered clustermap, for comparison with above
    hmap = sns.clustermap(
    df,
    method="average",
    metric="euclidean",
    col_cluster=False,
    row_cluster=False,
    yticklabels=False,
    cmap=cell_colourmap,
    col_colors=site_colours,
    row_colors=country_colours,
    cbar_kws={"label":"haplogroup"},
    cbar_pos=(1.05, 0.35, 0.05, 0.18)
    )
    add_clustermap_legends(hmap, len(cell_colourmap))
    hmap.savefig(f"{output_prefix}_hmap.pdf")

if __name__ == "__main__":
    main()
