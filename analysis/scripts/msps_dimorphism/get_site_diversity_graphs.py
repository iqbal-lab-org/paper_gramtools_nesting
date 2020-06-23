from typing import Union, List, NamedTuple, Optional
import sys
from pathlib import Path
import re
import json
from bisect import bisect_right

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
        - At least one site in the jvcf (by definition, that site must be a lvl1 site)
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
        print(f"No sites fall within specified region {region}")
        raise e

    next_idx = get_next_greater(cur_idx, lvl1_indices, num_sites)
    if next_idx == num_sites:
        next_idx = -1

    while True:
        print(cur_idx, next_idx)
        wire(cur_idx, next_idx, nesting_lvl, jvcf, result)
        cur_idx = next_idx
        next_idx = get_next_greater(cur_idx, lvl1_indices, num_sites)
        if not is_in_region(jvcf["Sites"][cur_idx], region):
            break
        if next_idx == num_sites:
            wire(cur_idx, -1, nesting_lvl, jvcf, result)
            break

    return result


@click.command()
@click.option(
    "--region",
    help="In the form 'SEG:start-end', as in samtools/bcftools",
    default=None,
    callback=get_region,
)
@click.argument("jvcf_input", type=click.Path(exists=True))
@click.argument("output_prefix", type=str)
def main(region: Region, jvcf_input: click.Path, output_prefix: str):
    output_dir = Path(output_prefix)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    else:
        assert output_dir.is_dir()

    with open(jvcf_input) as fin:
        jvcf = json.load(fin)


if __name__ == "__main__":
    main()
