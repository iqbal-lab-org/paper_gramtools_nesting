from typing import Union, List, Tuple, NamedTuple, Optional
import sys
from pathlib import Path
import re
import json
from bisect import bisect_right
from collections import defaultdict

import click
from igraph import Graph

region_matcher = re.compile(r"(\w+):(\d+)-(\d+)")


class Region(NamedTuple):
    segment: str = ""
    start: int = 0
    end: int = -1


def get_region(ctx, param, region_str: Optional[str]) -> Region:
    """
    callback validator to click command line parameter
    """
    if region_str is None:
        return Region()

    match_obj = region_matcher.fullmatch(region_str)
    if match_obj is not None:
        segment, start, end = (
            match_obj.group(1),
            int(match_obj.group(2)),
            int(match_obj.group(3)),
        )
        if start <= end and start >= 1 and end >= 1:
            return Region(segment, start, end)
    raise click.BadParameter(
        "region does not conform to requirements: 'seg:start-end', 1-based, end>=start."
    )


def get_sample_indices(jvcf, samples: List[str] = None) -> List[int]:
    result = list()
    # Map each sample name to the index it occurs at. This is also the index it occurs at in 'Sites'
    all_samples = {entry["Name"]: index for index, entry in enumerate(jvcf["Samples"])}
    if samples is None:
        result = list(range(len(all_samples)))
    else:
        for sample in samples:
            result.append(all_samples[sample])
    return result


def is_in_region(site_json, region: Region) -> bool:
    if region.segment == "":
        return True
    site_segment = site_json["SEG"]
    site_pos = int(site_json["POS"])
    if region.segment == site_segment:
        if region.start <= site_pos <= region.end:
            return True
    return False


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
        for site_indices in children.values():
            for array_idx, child_idx in enumerate(site_indices):
                if array_idx == 0:
                    result.add_edges([(cur_idx, child_idx)])
                if array_idx < len(site_indices) - 1:
                    wire(
                        child_idx,
                        site_indices[array_idx + 1],
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


def make_site_graph(jvcf, region: Region) -> Graph:
    """
    Assumptions:
        - There is at least one site in the jvcf (by definition, that site must be a lvl1 site)
    """
    result = Graph(directed=True)
    num_sites = len(jvcf["Sites"])
    nesting_lvl = 1

    # Reserve vertices
    result.add_vertices(num_sites)
    result.vs["POS"] = [0] * num_sites
    result.vs["populated"] = [False] * num_sites
    result.vs["nesting_lvl"] = [0] * num_sites

    lvl1_indices: List[int] = list()
    if jvcf["Lvl1_Sites"] == "all":
        lvl1_indices = list(range(num_sites))
    else:
        lvl1_indices = sorted(map(int, jvcf["Lvl1_Sites"]))

    # Position first index within :param: region
    cur_idx = 0
    try:
        while not is_in_region(jvcf["Sites"][cur_idx], region):
            cur_idx += 1
    except IndexError as e:
        print(f"ERROR: No sites fall within specified region {region}")
        raise IndexError from None

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
    result.vs["name"] = to_keep

    result["names"] = set(to_keep)

    return result


class DivReturn(NamedTuple):
    total_counts: int
    gt_heterozygosity: int


def compute_diversity(
    site_json, sample_indices: Optional[List[int]] = None
) -> DivReturn:
    num_samples = len(site_json["GT"])
    if sample_indices is None:
        sample_indices = list(range(num_samples))
    gt_counts = defaultdict(int)
    for idx in sample_indices:
        gts = site_json["GT"][idx]
        assert len(gts) == 1  # Assumes haploid for now
        if gts[0] is not None and "AMBIG" not in site_json["FT"]:
            gt_counts[gts[0]] += 1
    total_counts = sum(gt_counts.values())
    freqs = {gt: count / total_counts for gt, count in gt_counts.items()}

    diversity = 100 if total_counts > 0 else 0
    for freq in freqs.values():
        diversity -= freq * freq * 100

    return DivReturn(total_counts, int(diversity))


def annotate_vertices(
    graph: Graph,
    jvcf,
    sample_indices: List[int],
    attr_names: List[str],
    annotation_function,
):
    all_sites = jvcf["Sites"]
    new_attributes = {attr_name: list() for attr_name in attr_names}
    for idx, site_json in enumerate(all_sites):
        if idx not in graph["names"]:
            continue
        result = annotation_function(site_json, sample_indices)
        for attr_name in new_attributes:
            new_attributes[attr_name].append(getattr(result, attr_name))

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
@click.argument("output_dir", type=str)
def main(
    region: Region, partition_file: click.Path, jvcf_input: click.Path, output_dir: str
):
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    else:
        assert output_dir.is_dir()

    with open(jvcf_input) as fin:
        jvcf = json.load(fin)

    graph = make_site_graph(jvcf, region)

    if partition_file is None:
        annotate_vertices(
            graph, jvcf, None, ["total_counts", "gt_heterozygosity"], compute_diversity
        )
        graph.write_gml(str(output_dir / "all_samples.gml"))
    else:
        with open(partition_file) as fin:
            partitions = [line.rstrip().split("\t") for line in fin.readlines()]
            for idx, partition in enumerate(partitions):
                graph_copy = graph
                sample_indices = get_sample_indices(jvcf, partition)
                annotate_vertices(
                    graph_copy,
                    jvcf,
                    sample_indices,
                    ["total_counts", "gt_heterozygosity"],
                    compute_diversity,
                )
                graph.write_gml(str(output_dir / f"partition_{idx+1}.gml"))


if __name__ == "__main__":
    main()
