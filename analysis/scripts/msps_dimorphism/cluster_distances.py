import sys
from typing import Dict
from random import choice

import click
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode
from Bio import SeqIO
import edlib
from numpy import mean

from common import get_partition

NUM_SAMPLED_SEQS = 10000

SeqName = str
Seq = str
SeqMap = Dict[SeqName, Seq]


def get_seqmap(input_sequences: click.Path) -> SeqMap:
    result = dict()
    for record in SeqIO.parse(input_sequences, "fasta"):
        result[record.id] = str(record.seq)
    return result


def average_sampled_distance(cluster_1, cluster_2, seqmap: SeqMap):
    distances = list()
    for _ in range(NUM_SAMPLED_SEQS):
        seq1 = seqmap[choice(cluster_1)]
        seq2 = seqmap[choice(cluster_2)]
        distance = edlib.align(seq1, seq2, task="distance")["editDistance"] / len(seq1)
        distances.append(distance)
    return mean(distances)


@click.command()
@click.argument("hapg_matrix", type=click.Path(exists=True))
@click.argument("input_sequences", type=click.Path(exists=True))
@click.argument("output_prefix", type=str)
def main(hapg_matrix: click.Path, input_sequences: click.Path, output_prefix: str):
    clusters = get_partition(hapg_matrix)
    seqmap = get_seqmap(input_sequences)
    ofile = open(f"{output_prefix}_cluster_distmatrix.tsv", "w")
    for i in range(len(clusters)):
        row = list()
        for j in range(0, i):
            row.append("-")
        for j in range(i, len(clusters)):
            dist = average_sampled_distance(clusters[i], clusters[j], seqmap)
            row.append(str(dist))
        ofile.write("\t".join(row) + "\n")
    ofile.close()


if __name__ == "__main__":
    main()
