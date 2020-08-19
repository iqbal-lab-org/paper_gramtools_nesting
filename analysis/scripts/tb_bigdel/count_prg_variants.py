"""
"""

from pathlib import Path
from typing import List, Tuple
import json

import click
from pysam import VariantFile, VariantRecord

from tb_bigdel.common import Interval, Intervals, load_input_dels


class VarContainer(Interval):
    """
    An interval that holds variants
    """

    def __init__(self, start: int, stop: int, size: int):
        self.start = start
        self.stop = stop
        self.size = size
        self.num_vars, self.num_ref_bases = 0, 0

    @classmethod
    def make_from(cls, other: "Interval") -> "VarContainer":
        return VarContainer(other.start, other.stop, other.size)

    def __repr__(self):
        return f"[{self.start}, {self.stop}]"


VarContainers = List[VarContainer]


def add_if_spanned(interval: Interval, containers: VarContainers) -> None:
    """
    Record `interval` as spanned by a container
    """
    spanners = [cand for cand in containers if cand.spans(interval)]
    if len(spanners) > 1:
        raise ValueError(f"Found >1 variant spanning interval {interval}: {spanners}")
    elif len(spanners) == 1:
        if interval.spans(spanners[0]):
            # Case: the interval and the container are the same interval: we have found the input region and ignore it.
            return
        spanners[0].num_vars += 1
        spanners[0].num_ref_bases += len(interval)


def find_nested_ref_sites(json_prg, regions: VarContainers):
    """
    Because jvcf encodes nesting, need to only find sites in it
    that are nested within the ref section of each region in `regions`.

    I add a tolerance to the region searched because I have seen that about 5 deletion regions in gramtools vcf/jvcf have a start position a few bases after where the start position in the original bed.
    """
    TOLERANCE = 10
    result = []

    child_map = json_prg["Child_Map"]
    lvl1_sites = set(json_prg["Lvl1_Sites"])
    all_found = []

    for i, site in enumerate(json_prg["Sites"]):
        ref_allele = site["ALS"][0]
        site_interval = Interval(
            site["POS"] - TOLERANCE, site["POS"] + len(ref_allele) - 1
        )
        spanned = [reg for reg in regions if site_interval.spans(reg)]
        if len(spanned) == 0:
            continue
        elif len(spanned) > 1:
            raise ValueError(
                f"In jvcf, found interval {interval} spanning multiple input deletions"
            )

        # Remove any previously found interval that spans this found interval: this interval replaces it.
        all_found = [found for found in all_found if not found.spans(site_interval)]
        all_found.append(site_interval)
        # Get all sites directly under the found interval
        # Some might occur outside of the original region boundaries, they will get filtered out later
        nested_sites = [str(i)]
        added_sites = []
        while len(nested_sites) > 0:
            next_site = nested_sites.pop()
            for vals in child_map.get(next_site, dict()).values():
                nested_sites.extend(map(str, vals))
                for val in vals:
                    if str(val) not in child_map:
                        added_sites.append(val)
        print(
            f"Found sites {added_sites} under deletion site {site_interval} at idx {i}"
        )
        result += added_sites
    # Check for bijection between input regions and their representation in jvcf, warn if not
    missing_regions = []
    for reg in regions:
        if not any(map(lambda interval: interval.spans(reg), all_found)):
            missing_regions.append(reg)
    if len(missing_regions) > 0:
        print(f"Warning: could not find: {missing_regions} in jvcf")
    return result


def record_pair_yielder(vcf_file) -> Tuple[VariantRecord, VariantRecord]:
    vcf_recs = VariantFile(vcf_file)
    yielder = vcf_recs.fetch()
    cur_record = next(yielder)
    while True:
        try:
            next_record = next(yielder)
            yield (cur_record, next_record)
        except StopIteration:
            yield (cur_record, None)
            break
        cur_record = next_record


def add_vcf_metrics(vcf_file, input_dels: VarContainers) -> None:
    for cur_rec, next_rec in record_pair_yielder(vcf_file):
        cur_rec_stop = cur_rec.pos + len(cur_rec.ref) - 1
        # Ignore a record if it overlaps with next one
        if next_rec is not None and cur_rec_stop > next_rec.pos:
            print(cur_rec, next_rec)
            continue
        interval = Interval(cur_rec.pos, cur_rec_stop)
        add_if_spanned(interval, input_dels)


def add_jvcf_metrics(jvcf_file, input_dels: VarContainers) -> None:
    with open(jvcf_file) as fin:
        json_prg = json.load(fin)
    sites_to_consider = find_nested_ref_sites(json_prg, input_dels)
    for idx in sites_to_consider:
        json_site = json_prg["Sites"][idx]
        ref = json_site["ALS"][0]
        pos = json_site["POS"]
        interval = Interval(pos, pos + len(ref) - 1)
        add_if_spanned(interval, input_dels)


@click.command()
@click.argument(
    "variant_file", type=click.Path(exists=True),
)
@click.argument("input_dels_bed", type=click.Path(exists=True))
@click.argument("output_file", type=str)
def main(variant_file: click.Path, input_dels_bed: click.Path, output_file: str):
    """
    From an input bed of regions and a vcf/jvcf, count the number of variants and the number of ref variant bases under the regions.
    Gives a quantification of the diversity captured by a prg.

    :param: variant_file: a vcf or jvcf file describing variants a prg",
    """
    input_dels: Intervals = load_input_dels(input_dels_bed)
    input_dels = [VarContainer.make_from(reg) for reg in input_dels]

    suffixes = set(Path(variant_file).suffixes)
    if ".json" in suffixes:
        add_jvcf_metrics(variant_file, input_dels)
    elif len({".vcf", ".gz"} & suffixes) > 0:
        add_vcf_metrics(variant_file, input_dels)
    else:
        raise ValueError(f"{variant_file} has neither vcf nor jvcf file suffixes")
    with Path(output_file).open("w") as fout:
        print(f"region_start\tregion_len\tnum_vars\tbases_covered", file=fout)
        for input_del in input_dels:
            print(
                f"{input_del.start}\t{len(input_del)}\t{input_del.num_vars}\t{input_del.num_ref_bases}",
                file=fout,
            )


if __name__ == "__main__":
    main()
