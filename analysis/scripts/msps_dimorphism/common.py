from typing import List

import click
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode

country_to_colour = {"Ghana": "Gold", "Laos": "FireBrick", "Cambodia": "RoyalBlue"}
cluster_dimorphic_to_colour = {"form1": "Chocolate", "form2": "RosyBrown"}

def get_partition(hapg_matrix: click.Path) -> List[str]:
    df = pd.read_csv(hapg_matrix, sep="\t", index_col=0)
    cl = linkage(df, method="average", metric="euclidean")
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

