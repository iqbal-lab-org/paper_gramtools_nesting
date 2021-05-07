from typing import List
from collections import Counter, defaultdict

import click
import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode

country_to_colour = {"Ghana": "Gold", "Laos": "FireBrick", "Cambodia": "RoyalBlue"}
cluster_dimorphic_to_colour = {"form1": "Chocolate", "form2": "RosyBrown"}
nested_colour_mapping = {True: "gray", False: "navajowhite"}
marker_colours = {"Low": "Black", "High": "red"}

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

not_none = lambda x: x is not None

def get_complete_counts(gts1, gts2):
    """Express the genotype sets at the union of genotypes in the sets"""
    gt1_counts = Counter(filter(not_none,gts1))
    gt2_counts = Counter(filter(not_none,gts2))
    
    all_gts = set.union(set(gt1_counts.keys()),set(gt2_counts.keys()))
    gt1_distrib = [gt1_counts.get(gt,0) for gt in all_gts]
    gt2_distrib = [gt2_counts.get(gt,0) for gt in all_gts]
    if sum(gt1_distrib) == 0 and sum(gt2_distrib) == 0:
        raise ValueError(f"No non-null calls in {gts1} and {gts2}")
    return gt1_distrib, gt2_distrib

def allelic_distinguishability(gts1, gts2) -> float:
    """When randomly choosing two genotypes from the sets, how likely are they to be different?"""
    gt1_counts, gt2_counts = get_complete_counts(gts1, gts2)
    total_gt1 = sum(gt1_counts)
    total_gt2 = sum(gt2_counts)
    if total_gt1 == 0 or total_gt2 == 0:
        return 0
    else:
        prob_two_forms_same_alleles = 0
        for elem1, elem2 in zip(gt1_counts, gt2_counts):
            prob_two_forms_same_alleles += elem1 / total_gt1 * elem2 / total_gt2
        return 1 - prob_two_forms_same_alleles
    
def allelic_specificity(gts1, gts2) -> float:
    """How specific of a genotype set are the genotype calls?"""
    pre_none_ratio = len(gts1)/(len(gts1)+len(gts2))
    nonull_gts1 = list(filter(not_none,gts1))
    nonull_gts2 = list(filter(not_none,gts2))
    post_none_ratio = len(nonull_gts1)/(len(nonull_gts1) + len(nonull_gts2))
    return pre_none_ratio - post_none_ratio
