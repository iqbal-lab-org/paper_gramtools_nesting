"""
Converts a Bed file expressed in terms of a base reference into 
a Bed expressed in terms of a personalised reference.
"""
import sys
from csv import reader as csv_reader
from pathlib import Path

from gramtools.commands.genotype.seq_region_map import (
    SearchableSeqRegionsMap,
    BisectTarget,
    Chrom,
)


def translate_pos(chrom: Chrom, pos: int, searcher: SearchableSeqRegionsMap) -> int:
    region_idx = searcher.bisect(chrom, pos, BisectTarget.BASE_REF)
    region = searcher.get_region(chrom, region_idx)
    base_ref_offset = pos - region.base_ref_start
    translation = region.pers_ref_start + base_ref_offset
    return translation


def usage():
    print(f"{sys.argv[0]} input.bed rebasing_map.json output.bed")
    exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        usage()

    input_bed: Path = Path(sys.argv[1]).resolve()
    rebasing_map: Path = Path(sys.argv[2]).resolve()
    output_bed: Path = Path(sys.argv[3]).resolve()

    for fname in (input_bed, rebasing_map, output_bed.parent):
        if not fname.exists():
            print(f"Error: {fname} required but not found")
            usage()

    searcher = SearchableSeqRegionsMap.load_from(rebasing_map)
    with input_bed.open("r") as bed_in, output_bed.open("w") as bed_out:
        reader = csv_reader(bed_in, delimiter="\t")
        for line in reader:
            chrom = line[0]
            start_pos = (
                int(line[1]) - 1
            )  # Bed start is 0-based, SeqRegion coords are 1-based
            end_pos = int(line[2])

            translated_start_pos = translate_pos(chrom, start_pos, searcher)
            translated_start_pos = str(translated_start_pos - 1)
            translated_end_pos = translate_pos(chrom, end_pos, searcher).__str__()

            bed_out.write(
                "\t".join([chrom, translated_start_pos, translated_end_pos] + line[3:])
                + "\n"
            )
