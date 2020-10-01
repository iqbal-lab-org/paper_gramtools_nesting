"""
From a sam file and a genome containing the reference entries it refers to, find the
best alignment for each query sequence and calculate + output its Needleman-Wunsch edit
distance to the reference.
"""

import gzip
from typing import Dict, Tuple, List, Optional
from collections import defaultdict

from pysam import AlignmentFile, AlignedSegment
from Bio import SeqIO
import edlib
import click

Read = AlignedSegment
Chroms = Dict[str, str]

SAM_MATCH_OP = {0, 7, 8}
SAM_REF_CONSUMING_OP = {0, 2, 3, 7, 8}


def load_fasta(fname: str) -> Chroms:
    if fname.endswith("gz"):
        handle = gzip.open(fname, "rt")
    else:
        handle = open(fname, "r")

    result = dict()
    for record in SeqIO.parse(handle, "fasta"):
        result[record.id] = str(record.seq)
    handle.close()
    return result


class SamFileProcessingError(Exception):
    pass


def is_primary_alignment(read: Read) -> bool:
    """
    See SAM spec: 0x100 and 0x800 for secondary and supplementary alignment
    respectively, a single read with neither is primary alignment in sam file.
    """
    return read.flag & 0x900 == 0


def is_unmapped(read: AlignedSegment):
    return read.flag & 0x4 == 0x4


def is_match(cigar_tuple: Tuple) -> bool:
    return cigar_tuple[0] in SAM_MATCH_OP


def is_ref_consuming(cigar_tuple: Tuple) -> bool:
    return cigar_tuple[0] in SAM_REF_CONSUMING_OP


def get_ref_pos(read: AlignedSegment) -> int:
    result = read.reference_start  # pysam already converts it to 0-based
    for cig_t in read.cigartuples:
        if not is_ref_consuming(cig_t):
            result -= cig_t[1]
        else:
            break
    return result


def count_matches(read: AlignedSegment) -> int:
    result = 0
    if read.cigartuples is not None:
        for cig_t in read.cigartuples:
            if is_match(cig_t):
                result += cig_t[1]
    return result


def get_NW_edit_distance(query: str, target: str) -> int:
    alignment = edlib.align(query, target, mode="NW")
    return alignment["editDistance"]


class ReadMapping:
    def __init__(self):
        self.primary_alignment: Optional[Read] = None
        self.best_alignment: Optional[Read] = None
        self.num_matches: int = 0

    @property
    def seq(self):
        return self.primary_alignment.seq

    @property
    def name(self):
        return self.best_alignment.query_name

    @property
    def ref_name(self):
        return self.best_alignment.reference_name

    def update(self, other: Read):
        candidate_matches = count_matches(other)
        if self.best_alignment is None or candidate_matches > self.num_matches:
            self.best_alignment = other
            self.num_matches = candidate_matches
        if is_primary_alignment(other):
            if self.primary_alignment is not None:
                raise SamFileProcessingError(
                    f"Error: multiple primary alignments for query sequence {self.name}"
                )
            # Ensure the primary mapping contains the full query sequence, testing via
            # the CIGAR string whose inferred length is the full sequence
            candidate_length = len(other.seq)
            inferred_length = other.infer_read_length()
            if inferred_length is not None and candidate_length != inferred_length:
                raise SamFileProcessingError(
                    f"Primary alignment with sequence of length {candidate_length} has CIGAR string of length {inferred_length}"
                )
            self.primary_alignment = other

    def add_NM(self, chroms: Chroms):
        """
        Computes edit distance between self and the reference portion referred to by
        the alignment, and adds it as a tag to the read
        """
        if is_unmapped(self.best_alignment):
            return
        if self.ref_name not in chroms:
            raise ValueError(
                f"{self.ref_name} not in dictionary of chromosomes {chroms.keys()}"
            )
        reference_position = get_ref_pos(self.best_alignment)
        reference_portion = chroms[self.ref_name][
            reference_position : reference_position + len(self.seq)
        ]
        eddist = get_NW_edit_distance(self.seq, reference_portion)
        self.best_alignment.set_tag("NM", eddist)


BestMappings = Dict[str, ReadMapping]


@click.command()
@click.argument(
    "input_samfile",
    type=click.Path(exists=True),
)
@click.argument(
    "ref_genome_file",
    type=click.Path(exists=True),
)
@click.argument("output_samfile", type=str)
@click.option("--num_seqs", type=int, default=None)
def main(input_samfile, ref_genome_file, output_samfile, num_seqs):
    mappings: BestMappings = defaultdict(ReadMapping)
    sam_in = AlignmentFile(input_samfile, "r")
    for read in sam_in.fetch(until_eof=True):
        mappings[read.query_name].update(read)

    if num_seqs is not None and len(mappings) != num_seqs:
        raise SamFileProcessingError(
            f"Found {len(mappings)} distinct mappings but expected {num_seqs} from CLI"
        )
    chroms: Chroms = load_fasta(ref_genome_file)
    with AlignmentFile(output_samfile, "w", template=sam_in) as fout:
        for best_read in mappings.values():
            best_read.add_NM(chroms)
            fout.write(best_read.best_alignment)


if __name__ == "__main__":
    main()
