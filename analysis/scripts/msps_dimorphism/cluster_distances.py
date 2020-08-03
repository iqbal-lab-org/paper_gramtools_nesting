import sys
from typing import List, Dict
from random import choice

import click
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode
from Bio import SeqIO
import edlib
from numpy import mean

NUM_SAMPLED_SEQS = 10000


def get_partition(hapg_matrix: click.Path) -> List[str]:
    df = pd.read_csv(hapg_matrix, sep="\t", index_col=0)
    cl = linkage(df)
    root: ClusterNode = to_tree(cl)
    cur_node = root
    while True:
        cl_size = cur_node.get_count()
        child_sizes = [cur_node.left.get_count(), cur_node.right.get_count()]
        if sum(child_sizes) == 2:
            raise ValueError("No suitable partition of the cluster found")
        if min(child_sizes) < (cl_size / 10):
            largest = (
                cur_node.left if child_sizes[0] > child_sizes[1] else cur_node.right
            )
            cur_node = largest
        else:
            groups = [cur_node.left.pre_order(), cur_node.right.pre_order()]
            break
    sample_groups = []
    for group in groups:
        sample_names = [df.index[sample_idx] for sample_idx in group]
        sample_groups.append(sample_names)
    print(f"Found partition of sizes: {len(sample_groups[0])}, {len(sample_groups[1])}")
    assert len(set(sample_groups[0]).intersection(sample_groups[1])) == 0
    return sample_groups


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
