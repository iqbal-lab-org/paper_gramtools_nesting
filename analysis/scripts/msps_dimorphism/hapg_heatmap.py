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
    first_idx_in_region_non_nested,
)
from common import (
        get_partition, 
        country_to_colour, 
        nested_colour_mapping,
        marker_colours,
        allelic_distinguishability,
        allelic_specificity
)

SiteIdx = int
Hapgs = List[int]
Hapg_Dict = Dict[SiteIdx, Hapgs]

country_to_colour = {key: val.lower() for key, val in country_to_colour.items()}


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

    first_idx = first_idx_in_region_non_nested(jvcf, region)
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


def get_dimorphism_calls(hapg_matrix_file, region, jvcf):
    """ Partition genotype calls based on dimorphic form """
    groups=get_partition(hapg_matrix_file)
    sample_to_dimorphic = {sample_name: "form1" for sample_name in groups[0]}
    sample_to_dimorphic.update({sample_name: "form2" for sample_name in groups[1]})
    
    first_idx = first_idx_in_region_non_nested(jvcf, region)
        
    dimorphism_calls = defaultdict(list)
    cur_idx = first_idx
    while is_in_region(jvcf["Sites"][cur_idx], region):
        dimorphism_site_calls = defaultdict(list)
        site = jvcf["Sites"][cur_idx]
        gts = [gt[0] for gt in site["GT"]]
        for sample_idx, gt in enumerate(gts):
            sample_name = jvcf["Samples"][sample_idx]["Name"]
            try:
                form = sample_to_dimorphic[sample_name]
                dimorphism_site_calls[form].append(gt)
            except KeyError:
                pass

        for key,val in dimorphism_site_calls.items():
            dimorphism_calls[key].append(val)
        cur_idx += 1
    return dimorphism_calls

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

### Colourbars in clustermaps ###
def customise_cbar(hmap, num_colours):
    cbar = hmap.ax_heatmap.collections[0].colorbar
    r = cbar.vmax - cbar.vmin
    cbar.set_ticks(
        [
            cbar.vmin + 0.5 * r / (num_colours) + r * i / (num_colours)
            for i in range(num_colours)
        ]
    )
    cbar.set_ticklabels(list(range(-1, num_colours - 1)))
    
def add_nested_site_legend(hmap):
    for label,colour in nested_colour_mapping.items():
        hmap.ax_col_dendrogram.bar(0, 0, color=colour, label=label, linewidth=0)
    l1 = hmap.ax_col_dendrogram.legend(title='Nested site', loc="center", ncol=1, bbox_to_anchor=(0.6, 0.87), bbox_transform=gcf().transFigure)
    

def add_clustermap_legends(hmap, num_colours, countries = True):
    """See https://stackoverflow.com/q/27988846/12519542 for an excellent tutorial through this"""

    customise_cbar(hmap, num_colours)
    add_nested_site_legend(hmap)
    
    if countries:
        for label,colour in country_to_colour.items():
            hmap.ax_row_dendrogram.bar(0, 0, color=colour, label=label, linewidth=0)
        l2 = hmap.ax_row_dendrogram.legend(title='Sample country', loc="center", ncol=1, bbox_to_anchor=(0.25, 0.87), bbox_transform=gcf().transFigure)
    else:
        for label,colour in marker_colours.items():
            hmap.ax_row_dendrogram.bar(0, 0, color=colour, label=label, linewidth=0)
        l2 = hmap.ax_row_dendrogram.legend(title='Dimorphism specificity (top row)\nDimorphism sensitivity (bottom row)', loc="center left", 
                                           ncol=2, bbox_to_anchor=(0.05, 0.87), bbox_transform=gcf().transFigure)

    hmap.ax_heatmap.set(xlabel="Variant sites",ylabel="Samples")


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
    CBAR_POS = (1.04, 0.55, 0.05, 0.18)
    output_prefix = Path(output_prefix)

    with open(jvcf_input) as fin:
        jvcf = json.load(fin)


    hapgs_all_sites = get_hapgs_all_sites(jvcf, region)
    site_is_nested = hapgs_all_sites.pop("nested")
    sample_names = [sample["Name"] for sample in jvcf["Samples"]]
    df = pd.DataFrame(hapgs_all_sites, index=sample_names)
    hapg_matrix_file = f"{output_prefix}_hapgs.tsv" 
    df.to_csv(hapg_matrix_file, sep="\t")

    #### Get colourings ####
    ## Cell colouring: according to haplogroup
    cell_colourmap = get_cell_colouring(df)

    ## Row colouring: colour samples by country of origin
    country_colours = get_country_colouring(metadata_file, df)

    ## Column colouring: colour sites by whether they are nested
    site_is_nested = pd.Series(site_is_nested)
    nested_colours = site_is_nested.map(nested_colour_mapping)

    ## row-clustered clustermap
    hmap_clustered = sns.clustermap(
    df,
    method="average",
    metric="euclidean",
    col_cluster=False,
    yticklabels=False,
    cmap=cell_colourmap,
    col_colors=nested_colours,
    row_colors=country_colours,
    cbar_kws={"label":"haplogroup"},
    cbar_pos=CBAR_POS,
)
    add_clustermap_legends(hmap_clustered, len(cell_colourmap))
    hmap_clustered.savefig(f"{output_prefix}_hmap_clustered.pdf")
    hmap_clustered.savefig(f"{output_prefix}_hmap_clustered.svg")

    ## non-row-clustered clustermap, for comparison with above
    hmap = sns.clustermap(
    df,
    method="average",
    metric="euclidean",
    col_cluster=False,
    row_cluster=False,
    yticklabels=False,
    cmap=cell_colourmap,
    col_colors=nested_colours,
    row_colors=country_colours,
    cbar_kws={"label":"haplogroup"},
    cbar_pos=CBAR_POS,
    )
    add_clustermap_legends(hmap, len(cell_colourmap))
    hmap.savefig(f"{output_prefix}_hmap.pdf")

    ## clustermap with site-level dimorphism sensitivity and specificity ###
    dimorphism_calls = get_dimorphism_calls(hapg_matrix_file, region, jvcf)
    ## Compute measures of per-site dimorphism
    dimorphism_sensitivity = []
    dimorphism_specificity = []
    num_sites = len(list(dimorphism_calls.values())[0])
    for site_idx in range(num_sites):
        form1_gts = dimorphism_calls["form1"][site_idx]
        form2_gts = dimorphism_calls["form2"][site_idx]
        dimorphism_sensitivity.append(allelic_distinguishability(form1_gts,form2_gts))
        dimorphism_specificity.append(allelic_specificity(form1_gts,form2_gts))

    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    axs[0].hist(dimorphism_sensitivity)
    axs[0].set(title="Dimorphism sensitivity", xlabel="value",ylabel="Number of sites")
    axs[1].hist(dimorphism_specificity)
    axs[1].set(title="Dimorphism specificity", xlabel="value")
    fig.savefig(f"{output_prefix}_marker_histogram.pdf")

    dimo_sensi_colour = [marker_colours["High"] if val > 0.8 else marker_colours["Low"] for val in dimorphism_sensitivity]
    dimo_speci_colour = list()
    for val in dimorphism_specificity:
        if val > 0.3 or val < -0.3:
            dimo_speci_colour.append(marker_colours["High"])
        else:
            dimo_speci_colour.append(marker_colours["Low"])   

    hmap_with_markers = sns.clustermap(
        df,
        method="average",
        metric="euclidean",
        col_cluster=False,
        yticklabels=False,
        cmap=cell_colourmap,
        col_colors=[dimo_speci_colour, nested_colours, dimo_sensi_colour],
        cbar_kws={"label":"haplogroup"},
        cbar_pos=(1.04, 0.55, 0.05, 0.18)
    )
    add_clustermap_legends(hmap_with_markers, len(cell_colourmap), countries=False)
    hmap_with_markers.savefig(f"{output_prefix}_hmap_with_markers.pdf")

if __name__ == "__main__":
    main()
