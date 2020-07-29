from typing import Union, List, Tuple, NamedTuple, Optional, Dict
import re

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


def is_in_region(site_json, region: Region) -> bool:
    if region.segment == "":
        return True
    site_segment = site_json["SEG"]
    site_pos = int(site_json["POS"])
    if region.segment == site_segment:
        if region.start <= site_pos <= region.end:
            return True
    return False


def first_idx_in_region(sites_json, region: Region) -> int:
    result = 0
    try:
        while not is_in_region(sites_json[result], region):
            result += 1
    except IndexError as e:
        print(f"ERROR: No sites fall within specified region {region}")
        raise IndexError from None

    return result
