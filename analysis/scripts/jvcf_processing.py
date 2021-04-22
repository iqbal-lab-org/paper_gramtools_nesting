"""
Helper code to process JSON VCFs made by gramtools
"""
from typing import NamedTuple, Optional, Dict, Set, List, Union
import re
from collections import namedtuple

import click
import edlib

JVCF = Dict
SiteJson = Dict
SiteJsons = List[Dict]

### Define a region for filtering/selection purposes ###
region_matcher = re.compile(r"([^:]+):(\d+)-(\d+)")


class Region(NamedTuple):
    segment: str = ""
    start: int = 0
    end: int = -1


def get_region(region_str: Optional[str]) -> Region:
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
    raise ValueError(
        "region does not conform to requirements: 'seg:start-end', 1-based, end>=start."
    )


def click_get_region(ctx, param, region_str: Optional[str]) -> Region:
    """
    callback validator to click command line parameter;
    parses a string into a Region
    """
    try:
        return get_region(region_str)
    except ValueError as e:
        raise click.BadParameter(e)


def is_in_region(site_json: SiteJson, region: Region) -> bool:
    if region.segment == "":
        return True
    site_segment = site_json["SEG"]
    site_pos = int(site_json["POS"])
    if region.segment == site_segment:
        if region.start <= site_pos <= region.end:
            return True
    return False


def first_idx_in_region(sites_json: SiteJsons, region: Region) -> int:
    result = 0
    try:
        while not is_in_region(sites_json[result], region):
            result += 1
    except IndexError as e:
        print(f"ERROR: No sites fall within specified region {region}")
        raise IndexError from None

    return result


def get_sites_in_region(sites_json: SiteJsons, region: Region) -> SiteJsons:
    first_idx = first_idx_in_region(sites_json, region)
    last_idx = first_idx
    max_idx = len(sites_json)
    while last_idx < max_idx and is_in_region(sites_json[last_idx], region):
        last_idx += 1
    return sites_json[first_idx:last_idx]


def get_n_sites_starting_from_region(
    sites_json: SiteJsons, region: Region, num_sites: int
) -> SiteJsons:
    first_idx = first_idx_in_region(sites_json, region)
    total_num_sites = len(sites_json)
    if (first_idx + num_sites) > total_num_sites:
        raise ValueError(
            f"Getting {num_sites} starting from {region} requires more sites than available"
        )
    return sites_json[first_idx : first_idx + num_sites]


### Extracting information ###
def find_sample_index(jvcf: JVCF, sample_id: str) -> int:
    """
    Find which index in the Samples array a sample name occurs in
    """
    sample_index = -1
    for i in range(len(jvcf["Samples"])):
        if jvcf["Samples"][i]["Name"] == sample_id:
            sample_index = i
            break
    assert sample_index > -1
    return sample_index


class AlleleCall(NamedTuple):
    gt: Union[int, None]
    allele: str

    def has_call(self):
        return self.gt is not None


def get_indexed_attribute(site_json: SiteJson, attrib: str, index: int):
    return site_json[attrib][index]


def get_called_allele(sample_index: int, site_json: SiteJson) -> AlleleCall:
    """
    Take the genotype call (GT) of a sample and extract the called allele from site's list of alleles (ALS)
    """
    allele_index = get_indexed_attribute(site_json, "GT", sample_index)
    assert len(allele_index) > 0
    allele_index = allele_index[0]

    allele = ""
    if allele_index is not None:
        allele = get_indexed_attribute(site_json, "ALS", allele_index)
    return AlleleCall(allele_index, allele)


def is_nested(lvl1sites: Set, site_idx: int) -> bool:
    return site_idx not in lvl1sites


def num_sites_under(child_map: Dict, site_idx: str) -> int:
    result = 0
    to_visit: List[int] = [site_idx]
    while len(to_visit) > 0:
        child_entries = child_map.get(to_visit.pop(), dict())
        for vals in child_entries.values():
            result += len(vals)
            to_visit.extend(map(str, vals))
    return result


### Convert a jVCF to a VCF ###
class jVCF_to_VCF:
    def __init__(self):
        pass

    def convert(self, jvcf: JVCF, fout):
        headers = self.make_headers(jvcf)
        fout.write(headers)
        for site in jvcf["Sites"]:
            fout.write(self.convert_one_site(site))

    def make_headers(self, jvcf) -> str:
        default_double_headers = (
            "##fileformat=VCFv4.2\n"
            '##FILTER=<ID=PASS,Description="All filters passed">\n'
            "##source=jVCF\n"
        )
        site_field_headers = ""
        for sf_name, sf in jvcf["Site_Fields"].items():
            tagname = "FORMAT" if sf_name != "FT" else "FILTER"
            site_field_headers += (
                f"##{tagname}=<ID={sf_name},Description=\"{sf['Desc']}\">\n"
            )

        samples = "\t".join([sample["Name"] for sample in jvcf["Samples"]])
        sample_header = (
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
        )
        return f"{default_double_headers}{site_field_headers}{sample_header}"

    def get_sample_fields(self, site_json: SiteJson) -> List:
        site_wide_fields = {"SEG", "POS", "ALS", "DP"}
        sample_fields = sorted(set.difference(set(site_json.keys()), site_wide_fields))
        # Place 'GT' at first position of FORMAT, as per VCF convention
        gt_loc = 0
        for i, field in enumerate(sample_fields):
            if field == "GT":
                gt_loc = i
                break
        if gt_loc != 0:
            tmp = sample_fields[0]
            sample_fields[0] = "GT"
            sample_fields[gt_loc] = tmp
        if len(sample_fields) == 0:
            raise ValueError(
                f"No sample-specific fields (not {site_wide_fields}) in site {site_json}"
            )
        return sample_fields

    def convert_one_site(self, site_json: SiteJson) -> str:
        result = f'{site_json["SEG"]}\t{site_json["POS"]}\t.\t'
        alleles = site_json["ALS"]
        if len(alleles) == 0:
            raise ValueError(f"No alleles in site {site_json}")
        result += f"{alleles[0]}\t"
        if len(alleles) == 1:
            result += ".\t"
        else:
            result += ",".join(alleles[1:]) + "\t"
        result += ".\t.\t.\t"  # QUAL, FILTER, INFO, currently not used

        sample_fields = self.get_sample_fields(site_json)
        result += ":".join(sample_fields) + "\t"  # FORMAT

        num_samples = len(site_json[sample_fields[0]])
        sample_entries = [[] for _ in range(num_samples)]
        for sample_field in sample_fields:
            next_field = site_json[sample_field]
            for i, sample_value in enumerate(next_field):
                used_value = sample_value
                if (
                    type(sample_value) is not list
                ):  # Some fields have multiple entries per sample, others not; make it uniform.
                    used_value = [used_value]
                if sample_field == "FT" and len(used_value) == 0:
                    used_value = ["PASS"]
                elif sample_field == "COV":
                    used_value = [int(val) for val in used_value]
                if sample_field == "GT":
                    used_value = ["." if val is None else val for val in used_value]
                    sample_entries[i].append("/".join(map(str, used_value)))
                else:
                    sample_entries[i].append(",".join(map(str, used_value)))
        serialised_sample_entries = "\t".join(
            map(lambda entry: ":".join(entry), sample_entries)
        )
        result += f"{serialised_sample_entries}\n"
        return result


### Evaluate a genotyped jvcf using a truth jvcf ###
class FixedDict(dict):
    """Cannot add new keys"""

    def __init__(self, other):
        dict.__init__(self)
        for key in other.keys():
            dict.__setitem__(self, key, other[key])

    def __setitem__(self, key, item):
        if key not in self:
            raise KeyError(f"Key {key} is not defined")
        dict.__setitem__(self, key, item)


def call_classif(res_call: AlleleCall, truth_call: AlleleCall) -> str:
    if res_call.has_call():
        if not truth_call.has_call():
            return "FP"
        else:
            if res_call.allele == truth_call.allele:
                return "TP"
            else:
                return "FP"
    else:
        if truth_call.has_call():
            return "FN"
        else:
            return "TN"


def evaluate_site(
    genotyped: SiteJson,
    genotyped_sample_index: int,
    truth: SiteJson,
    truth_sample_index: int,
) -> Dict:
    result_fields = [
        "res_has_call",
        "truth_has_call",
        "res_is_correct",
        "classif",
        "GC",
        "GCP",
        "edit_dist",
        "truth_allele",
        "gt_allele",
        "cov_gt_allele",
        "cov_other_alleles",
        "genotyped_ambiguous",
        "truth_ambiguous",
    ]
    result = FixedDict({k: "NA" for k in result_fields})
    result["genotyped_ambiguous"] = 0
    result["truth_ambiguous"] = 0

    truth_call = get_called_allele(truth_sample_index, truth)

    res_call = get_called_allele(genotyped_sample_index, genotyped)
    result["res_has_call"] = res_call.has_call()
    result["truth_has_call"] = truth_call.has_call()
    result["res_is_correct"] = res_call.allele == truth_call.allele
    result["classif"] = call_classif(res_call, truth_call)

    result["GC"] = get_indexed_attribute(genotyped, "GT_CONF", genotyped_sample_index)
    result["GCP"] = get_indexed_attribute(
        genotyped, "GT_CONF_PERCENTILE", genotyped_sample_index
    )

    result["truth_allele"] = truth_call.allele
    result["gt_allele"] = res_call.allele
    result["edit_dist"] = edlib.align(res_call.allele, truth_call.allele)[
        "editDistance"
    ]

    total_cov = get_indexed_attribute(genotyped, "DP", genotyped_sample_index)
    result["cov_other_alleles"] = total_cov
    if res_call.has_call():
        cov_gt_allele = get_indexed_attribute(genotyped, "COV", genotyped_sample_index)[
            res_call.gt
        ]
        result["cov_gt_allele"] = cov_gt_allele
        result["cov_other_alleles"] -= cov_gt_allele

    if "AMBIG" in get_indexed_attribute(genotyped, "FT", genotyped_sample_index):
        result["genotyped_ambiguous"] = 1
    if "AMBIG" in get_indexed_attribute(truth, "FT", truth_sample_index):
        result["truth_ambiguous"] = 1

    return result
