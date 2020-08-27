from pathlib import Path
from typing import Tuple
import sys

from pysam import AlignmentFile


def usage():
    print(
        f"usage: {sys.argv[0]} input_dir output_file\n"
        "The input_dir should contain the .sam files to analyse.\n"
    )
    exit(1)


def get_sample_and_condition_name(sam_fname: Path) -> Tuple[str, str]:
    elements = sam_fname.stem.split("_")
    sample = elements[-1]
    condition = "_".join(elements[:-1])
    return sample, condition


if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()

    input_dir = Path(sys.argv[1]).resolve()
    if not input_dir.exists():
        print(f"Error: {input_dir} not found")
        usage()

    output_bed = Path(sys.argv[2]).resolve()
    output_bed.parent.mkdir(exist_ok=True)

    sam_file_list = list(input_dir.glob(f"*.sam"))
    if len(sam_file_list) == 0:
        print(f"Error: no .sam files in {input_dir}")
        usage()

    with output_bed.open("w") as fout:
        for sam_fname in sam_file_list:
            sample, condition = get_sample_and_condition_name(sam_fname)
            samfile = AlignmentFile(sam_fname, "r")
            for read in samfile.fetch(until_eof=True):
                start = int(read.pos) - 1
                if start < 0:
                    continue
                stop = start + len(read.seq)
                qname = read.qname
                fout.write(f"{sample}\t{start}\t{stop}\t{condition}\t{qname}\n")
