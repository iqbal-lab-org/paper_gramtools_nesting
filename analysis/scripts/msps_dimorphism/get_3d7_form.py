import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode, dendrogram
from common import get_partition, country_to_colour, cluster_dimorphic_to_colour
from Bio import SeqIO
import edlib
import click

import matplotlib.pyplot as plt

def get_seq(fname):
    for seq_record in SeqIO.parse(fname, "fasta"):
        return str(seq_record.seq)

@click.command()
@click.argument("3d7_seq_fname",type=click.Path(exists=True))
@click.argument("sample_seq_dirname",type=click.Path(exists=True))
@click.argument("hapg_matrix_fname",type=click.Path(exists=True))
@click.argument("output_dirname",type=click.Path(exists=True))
def main(
        3d7_seq_fname,
        sample_seq_dirname,
        hapg_matrix_fname,
        output_dirname
):
    ref = get_seq(3d7_seq_fname)
    #ref = get_seq(f"../../../tmp_work/heatmaps/3D7_DBL_DBLMSP2.fa")
    #hapg_matrix="../../../tmp_work/heatmaps/DBL_DBLMSP2_hapgs.tsv"

    groups=get_partition(hapg_matrix_fname)

    distances = dict()
    group1=set(groups[0])
    group2=set(groups[1])
    for sname in groups[0] + groups[1]:
        seq = get_seq(f"{sample_seq_dirname}/{sname}.fa")
        distance = edlib.align(ref, seq)["editDistance"]
        distances[sname] = distance

    form1_dists = [distances[sname] for sname in distances if sname in group1]
    form2_dists = [distances[sname] for sname in distances if sname in group2]
    assert len(form1_dists + form2_dists) == len(groups[0] + groups[1])

    min_form1 = min(form1_dists)
    min_form2 = min(form2_dists)
    print(f'Min 3D7 distance to form1 samples: {min_form1}')
    print(f'Min 3D7 distance to form2 samples: {min_form2}')

    df = pd.read_csv(hapg_matrix_fname, sep="\t", index_col=0)
    cl = linkage(df, method="average", metric="euclidean")
    dn = dendrogram(cl)

    ## The dendrogram stores leaf labels as the index in the original data (df),
    ## left-to-right along it.
    ## Below goes from leaf label to sample name to distance, and if the distance is the closest to 3d7
    ## of the whole set, we give their corresponding dendrogram leaf a label saying so
    global_min = min(min_form1, min_form2)
    belongings = []
    original_positions = dict()
    for i, sname in enumerate(df.index):
        original_positions[sname] = i
    leaf_labels = ["" for _ in range(len(df.index))]
    for idx in dn["ivl"]: # traverses dendrogram leaves
        sname = df.index[int(idx)]
        assert sname in group1 or sname in group2
        if sname in group1:
            belongings.append("group1")
        else:
            belongings.append("group2")
        if distances[sname] == global_min:
            original_index = original_positions[sname]
            leaf_labels[original_index] = "REF_CLOSEST"

    # print() confirms that the groups extracted above are contiguous in the dendrogram
    #print(belongings)

    # Note: dendrogram must look the same as that in seaborn heatmap, confirm this visually
    fig = plt.figure(figsize=(25, 12))
    dn = dendrogram(cl,labels=leaf_labels)
    plt.savefig("{output_dirname}/dendrogram_ref_closest.pdf")
