from typing import Union, List, Optional, Dict
import sys
from pathlib import Path
import json
from math import log
from bisect import bisect_right
from collections import defaultdict
from itertools import combinations

import click
from igraph import Graph

from jvcf_processing import (
    Region,
    get_region,
    is_in_region,
    first_idx_in_region,
)


Sample_Indices = List[int]
FreqDict = Dict[str, float]
Distrib = List[float]


def get_sample_indices(jvcf, samples: List[str] = None) -> Sample_Indices:
    result = list()
    # Map each sample name to the index it occurs at. This is also the index it occurs at in 'Sites'
    all_samples = {entry["Name"]: index for index, entry in enumerate(jvcf["Samples"])}
    if samples is None:
        result = list(range(len(all_samples)))
    else:
        for sample in samples:
            result.append(all_samples[sample])
    return result


def wire(cur_idx: int, next_idx: int, nesting_lvl: int, jvcf, result: Graph) -> None:
    """
    Add vertices and edges to the :param: result graph using the nesting structure of :param: jvcf.
    Is a recursive function, calling itself to wire child(=nested) sites whenever they are encountered.
    """
    if not result.vs["populated"][cur_idx]:
        result.vs[cur_idx]["POS"] = jvcf["Sites"][cur_idx]["POS"]
        result.vs[cur_idx]["nesting_lvl"] = nesting_lvl
        result.vs[cur_idx]["populated"] = True

    if next_idx <= 0:
        return

    child_map = jvcf["Child_Map"]
    if str(cur_idx) not in child_map:
        result.add_edges([(cur_idx, next_idx)])
    else:
        children = jvcf["Child_Map"][str(cur_idx)]
        for child_indices in children.values():
            # it is vital to sort the child indices, as they are numbered according to their POS (lowest first), and we want graph topology to reflect POS
            sorted_child_indices = sorted(child_indices)
            for array_idx, child_idx in enumerate(sorted_child_indices):
                if array_idx == 0:
                    result.add_edges([(cur_idx, child_idx)])
                if array_idx < len(child_indices) - 1:
                    wire(
                        child_idx,
                        sorted_child_indices[array_idx + 1],
                        nesting_lvl + 1,
                        jvcf,
                        result,
                    )
                else:
                    wire(child_idx, next_idx, nesting_lvl + 1, jvcf, result)


def get_next_greater(idx: int, idx_list: List[int], max_val: int) -> int:
    insertion_point = bisect_right(idx_list, idx)
    if insertion_point == len(idx_list):
        return max_val
    else:
        return idx_list[insertion_point]


def new_site_graph(num_sites) -> Graph:
    result = Graph(directed=True)
    # Reserve vertices
    result.add_vertices(num_sites)
    result.vs["POS"] = [0] * num_sites
    result.vs["populated"] = [False] * num_sites
    result.vs["nesting_lvl"] = [0] * num_sites

    return result


def site_graph_copy(other: Graph) -> Graph:
    result = new_site_graph(len(other.vs))
    for attr in result.vs[0].attributes():
        result.vs[attr] = other.vs[attr]
    for edge in other.es:
        result.add_edges([(edge.source, edge.target)])
    result["idxs_in_prg"] = other["idxs_in_prg"]
    result.vs["idx_in_prg"] = other.vs[
        "idx_in_prg"
    ]  # To keep info on absolute site index
    return result


def make_site_graph(jvcf, region: Region) -> Graph:
    """
    Assumptions:
        - There is at least one site in the jvcf (by definition, that site must be a lvl1 site)
    """
    num_sites = len(jvcf["Sites"])
    nesting_lvl = 1

    result = new_site_graph(num_sites)

    lvl1_indices: List[int] = list()
    if jvcf["Lvl1_Sites"] == "all":
        lvl1_indices = list(range(num_sites))
    else:
        lvl1_indices = sorted(map(int, jvcf["Lvl1_Sites"]))

    cur_idx = first_idx_in_region(jvcf["Sites"], region)

    next_idx = get_next_greater(cur_idx, lvl1_indices, num_sites)
    if next_idx == num_sites:
        next_idx = -1

    while True:
        wire(cur_idx, next_idx, nesting_lvl, jvcf, result)
        cur_idx = next_idx
        next_idx = get_next_greater(cur_idx, lvl1_indices, num_sites)
        if next_idx == num_sites:
            # Mark the last site as processed
            wire(cur_idx, -1, nesting_lvl, jvcf, result)
            break
        if not is_in_region(jvcf["Sites"][cur_idx], region):
            break

    # Clear out non-populated vertices
    to_delete = list()
    to_keep = list()
    for idx, node in enumerate(result.vs):
        if not node["populated"]:
            to_delete.append(idx)
        else:
            to_keep.append(idx)
    result.delete_vertices(to_delete)
    result.vs["idx_in_prg"] = to_keep

    result["idxs_in_prg"] = set(to_keep)

    return result


diversity_metrics = {
    "gt_non_null_counts": 0,
    "gt_heterozygosity": 0.0,
    "hapg_heterozygosity": 0.0,
    "num_ambig": 0,
}


def heterozygosity(freqs: Dict[str, float]) -> float:
    """Computes heterozygosity, prob that two randomly sampled individuals differ at a site"""

    diversity = 1.0 if len(freqs) > 0 else 0.0
    for freq in freqs.values():
        diversity -= freq * freq
    if diversity == 0:
        diversity = 0.001  # To force use of float by cytoscape
    return diversity


def counts_to_freqs(counts_dict) -> Dict:
    total_counts = sum(counts_dict.values())
    freqs = {key: count / total_counts for key, count in counts_dict.items()}
    return freqs


def get_site_freqs(site_json, sample_indices: Sample_Indices):
    """Computes hapg and gt frequencies in a site, and also number of AMBIG sites, for a given set of samples"""
    gt_counts, hapg_counts = defaultdict(int), defaultdict(int)
    num_ambigs = 0
    for idx in sample_indices:
        gts, hapgs = site_json["GT"][idx], site_json["HAPG"][idx]
        assert len(gts) == 1  # Assumes haploid for now
        if "AMBIG" in site_json["FT"][idx]:
            num_ambigs += 1
            continue
        if gts[0] is not None:
            gt_counts[gts[0]] += 1
            if len(hapgs) > 0 and hapgs[0] is not None:
                hapg_counts[hapgs[0]] += 1
    return (
        counts_to_freqs(gt_counts),
        counts_to_freqs(hapg_counts),
        num_ambigs,
        sum(gt_counts.values()),
    )


def compute_diversity(site_json, *partitions: List[Optional[Sample_Indices]]) -> Dict:
    assert len(partitions) == 1
    result = diversity_metrics.copy()
    num_samples = len(site_json["GT"])
    sample_indices = partitions[0]
    if sample_indices is None:
        sample_indices = list(range(num_samples))

    gt_freqs, hapg_freqs, num_ambigs, counts = get_site_freqs(site_json, sample_indices)
    result["gt_heterozygosity"] = heterozygosity(gt_freqs)
    result["hapg_heterozygosity"] = heterozygosity(hapg_freqs)
    result["num_ambig"] = num_ambigs
    result["gt_non_null_counts"] = counts

    return result


divergence_metrics = {
    "gtJensen_Shannon": 0,
    "hapgJensen_Shannon": 0,
    "p1_gt_non_null_counts": 0,
    "p2_gt_non_null_counts": 0,
    "p1_ambig_counts": 0,
    "p2_ambig_counts": 0,
    "p1_hapg_distrib": "",
    "p2_hapg_distrib": "",
}


def normalise_distrib(freqs: FreqDict, all_keys: List[str]) -> Distrib:
    """Converts :param: freqs into a distribution with one entry per entry in :param: all_keys"""
    result = [0] * len(all_keys)
    for idx, key in enumerate(all_keys):
        result[idx] = freqs.get(key, 0)
    return result


def Kullback_Leibler(d1: Distrib, d2: Distrib) -> float:
    assert len(d1) == len(d2)
    result = 0
    for v1, v2 in zip(d1, d2):
        if v1 != 0 and v2 != 0:
            result += v1 * log(v1 / v2, 2)
    return result


def Jensen_Shannon(freqs_1: FreqDict, freqs_2: FreqDict) -> float:
    all_keys = set(freqs_1.keys()).union(set(freqs_2.keys()))
    all_keys = sorted(list(all_keys))  # Make sure the gt/hapg are output in order
    normed_1 = normalise_distrib(freqs_1, all_keys)
    normed_2 = normalise_distrib(freqs_2, all_keys)

    middle_distrib = [(normed_1[i] + normed_2[i]) / 2 for i in range(len(all_keys))]
    jensen_shannon = (
        Kullback_Leibler(normed_1, middle_distrib)
        + Kullback_Leibler(normed_2, middle_distrib)
    ) / 2
    if jensen_shannon == 0:
        jensen_shannon = 0.001  # Force float for cytoscape
    return (jensen_shannon, normed_1, normed_2)


def compute_Jensen_Shannon(site_json, *partitions: List[Sample_Indices]) -> Dict:
    assert len(partitions) == 2
    result = divergence_metrics.copy()
    gtfreqs_1, hapgfreqs_1, ambig_counts1, gt_counts1 = get_site_freqs(
        site_json, partitions[0]
    )
    gtfreqs_2, hapgfreqs_2, ambig_counts2, gt_counts2 = get_site_freqs(
        site_json, partitions[1]
    )

    result["gtJensen_Shannon"], _, _ = Jensen_Shannon(gtfreqs_1, gtfreqs_2)
    (result["hapgJensen_Shannon"], normed_1, normed_2,) = Jensen_Shannon(
        hapgfreqs_1, hapgfreqs_2
    )

    result["p1_hapg_distrib"] = str(normed_1)
    result["p2_hapg_distrib"] = str(normed_2)

    result["p1_gt_non_null_counts"] = gt_counts1
    result["p2_gt_non_null_counts"] = gt_counts2

    result["p1_ambig_counts"] = ambig_counts1
    result["p2_ambig_counts"] = ambig_counts2

    return result


def annotate_vertices(
    graph: Graph,
    jvcf,
    partitions: List[Sample_Indices],
    attr_names: List[str],
    annotation_function,
):
    all_sites = jvcf["Sites"]
    new_attributes = {attr_name: list() for attr_name in attr_names}
    for idx, site_json in enumerate(all_sites):
        if idx not in graph["idxs_in_prg"]:
            continue
        result = annotation_function(site_json, *partitions)
        for attr_name, value in result.items():
            new_attributes[attr_name].append(value)

    for attr_name, attr_values in new_attributes.items():
        graph.vs[attr_name] = attr_values


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

    graph = make_site_graph(jvcf, region)

    if partition_file is None:
        graph_copy = site_graph_copy(graph)
        annotate_vertices(
            graph_copy, jvcf, [None], list(diversity_metrics.keys()), compute_diversity
        )
        graph_copy.write_gml(f"{output_prefix}_all_samples.gml")
    else:
        with open(partition_file) as fin:
            partitions = [line.rstrip().split("\t") for line in fin.readlines()]
            # Diversity metrics for each partition
            for idx, partition in enumerate(partitions):
                graph_copy = site_graph_copy(graph)
                sample_indices = get_sample_indices(jvcf, partition)
                annotate_vertices(
                    graph_copy,
                    jvcf,
                    [sample_indices],
                    list(diversity_metrics.keys()),
                    compute_diversity,
                )
                graph_copy.write_gml(f"{output_prefix}_partition_{idx+1}.gml")

            # Divergence metrics for each pair of partitions
            for combination in combinations(range(len(partitions)), 2):
                partition1 = get_sample_indices(jvcf, partitions[combination[0]])
                partition2 = get_sample_indices(jvcf, partitions[combination[1]])

                graph_copy = site_graph_copy(graph)
                annotate_vertices(
                    graph_copy,
                    jvcf,
                    [partition1, partition2],
                    list(divergence_metrics.keys()),
                    compute_Jensen_Shannon,
                )

                combin_name = f"partition_{combination[0]+1}_vs_{combination[1]+1}.gml"
                graph_copy.write_gml(f"{output_prefix}_{combin_name}")


if __name__ == "__main__":
    main()
