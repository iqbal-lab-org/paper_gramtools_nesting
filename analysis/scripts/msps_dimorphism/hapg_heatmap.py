from typing import Dict, List
from pathlib import Path
from collections import defaultdict
import json

import click
import pandas as pd
import seaborn as sns

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
    is 0 by convention.
    """
    result: Hapgs = [0] * num_samples
    for sample_idx, gt in enumerate(site_json["GT"]):
        if gt[0] is not None:
            sample_hapg = site_json["HAPG"][sample_idx][0]
            result[sample_idx] = sample_hapg + 1
    return result


def get_hapgs_all_sites(jvcf, region: Region) -> Hapg_Dict:
    num_sites = len(jvcf["Sites"])
    num_samples = len(jvcf["Samples"])
    lvl1_sites = set(jvcf["Lvl1_Sites"])
    result: Hapg_Dict = dict()

    first_idx = first_idx_in_region(jvcf["Sites"], region)
    # If the first idx is not at lvl1, we can get a site following the first one which has a smaller POS and that is not in target region
    while first_idx not in lvl1_sites:
        first_idx += 1
    cur_idx = first_idx

    while is_in_region(jvcf["Sites"][cur_idx], region):
        result[cur_idx - first_idx] = get_hapgs_one_site(
            jvcf["Sites"][cur_idx], num_samples
        )
        cur_idx += 1
        if cur_idx == num_sites:
            break

    return result


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
    df = pd.DataFrame(hapgs_all_sites)
    df.to_csv(f"{output_prefix}_hapgs.tsv", sep="\t", index=False)

    hmap = sns.clustermap(df, col_cluster=False)
    hmap.savefig(f"{output_prefix}_hmap.pdf")


if __name__ == "__main__":
    main()
