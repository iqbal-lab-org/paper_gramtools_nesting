import sys
from pathlib import Path
from typing import List, NamedTuple

from ete3 import Tree, TreeStyle, RectFace, NodeStyle, TreeNode, add_face_to_node
import click
import pandas as pd

from common import get_partition, country_to_colour, cluster_dimorphic_to_colour

sample_to_country = dict()
phylo_colouring = dict()
cluster_colouring = dict()


def get_cluster_dimorphic_colours(groups):
    cluster_sample_to_dimorphic = {sample_name: "form1" for sample_name in groups[0]}
    cluster_sample_to_dimorphic.update(
        {sample_name: "form2" for sample_name in groups[1]}
    )
    result = dict()
    for sname, form in cluster_sample_to_dimorphic.items():
        result[sname] = cluster_dimorphic_to_colour[form]
    return result


def get_phylo_dimorphic_colours(tree: Tree, groups):
    cur_node = tree
    while True:
        branching = [node for node in cur_node.children if not node.is_leaf()]
        if len(branching) == 2:
            break
        elif len(branching) == 0:
            raise ValueError("Could not find phylo tree bifurcation")
        else:
            cur_node = branching[0]

    phylo_sample_to_dimorphic = {
        sample_name: "form1" for sample_name in branching[0].get_leaf_names()
    }
    phylo_sample_to_dimorphic.update(
        {sample_name: "form2" for sample_name in branching[1].get_leaf_names()}
    )

    ## Give same colour to majority form across the methods
    same_form_size = len(
        set(branching[0].get_leaf_names()).intersection(set(groups[0]))
    )
    diff_form_size = len(
        set(branching[0].get_leaf_names()).intersection(set(groups[1]))
    )
    if diff_form_size > same_form_size:
        phylo_dimorphic_to_colour = {
            "form1": cluster_dimorphic_to_colour["form2"],
            "form2": cluster_dimorphic_to_colour["form1"],
        }
    else:
        phylo_dimorphic_to_colour = cluster_dimorphic_to_colour.copy()

    result = dict()
    for sname, form in phylo_sample_to_dimorphic.items():
        result[sname] = phylo_dimorphic_to_colour[form]
    return result


def node_layout(node):
    if node.is_leaf():
        country_colour = country_to_colour[sample_to_country[node.name]]
        country_colour_face = RectFace(
            width=40, height=15, fgcolor=None, bgcolor=country_colour
        )
        add_face_to_node(country_colour_face, node, column=0, aligned=True)

        try:
            dimorphic_colour = phylo_colouring[node.name]
            dimorphic_colour_face = RectFace(
                width=40, height=15, fgcolor=None, bgcolor=dimorphic_colour
            )
            dimorphic_colour_face.margin_left = 10
            add_face_to_node(dimorphic_colour_face, node, column=1, aligned=True)
        except KeyError:
            pass

        try:
            dimorphic_colour = cluster_colouring[node.name]
            dimorphic_colour_face = RectFace(
                width=40, height=15, fgcolor=None, bgcolor=dimorphic_colour
            )
            dimorphic_colour_face.margin_left = 10
            add_face_to_node(dimorphic_colour_face, node, column=2, aligned=True)
        except KeyError:
            pass

    node_style = NodeStyle()
    node_style["size"] = 0
    node_style["vt_line_width"] = node_style["vt_line_width"] = 2
    node.set_style(node_style)


def get_treestyle():
    """This needs to be called multiple times to .render() multiple trees"""
    ts = TreeStyle()
    # ts.scale=120
    ts.show_leaf_name = False
    ts.arc_start = 0
    ts.arc_span = 360
    ts.mode = "c"
    ts.layout_fn = node_layout
    return ts


@click.command()
@click.argument("hapg_matrix_file", type=click.Path(exists=True))
@click.argument("tree_file", type=click.Path(exists=True))
@click.argument("sample_metadata_file", type=click.Path(exists=True))
@click.argument("output_prefix", type=str)
def main(
    hapg_matrix_file: click.Path,
    tree_file: click.Path,
    sample_metadata_file: click.Path,
    output_prefix: str,
):
    t = Tree(tree_file)

    sample_names = t.get_leaf_names()

    groups = get_partition(hapg_matrix_file)
    # Below make sure we have the same set of samples; only works if partition includes all samples
    # assert set(groups[0]).union(set(groups[1])) == set(sample_names)

    global cluster_colouring
    cluster_colouring = get_cluster_dimorphic_colours(groups)

    global phylo_colouring
    phylo_colouring = get_phylo_dimorphic_colours(t, groups)

    metadata = pd.read_csv(sample_metadata_file, sep="\t")
    global sample_to_country
    sample_to_country = dict(zip(list(metadata["sample"]), list(metadata["country"])))

    base_path = str(Path(output_prefix).resolve())
    t.render(f"{base_path}_tree.pdf", w=120, units="mm", tree_style=get_treestyle())

    ts = get_treestyle()
    ts.force_topology = True
    t.render(
        f"{base_path}_tree_fixed_branch_lengths.pdf", w=120, units="mm", tree_style=ts
    )


def draw_tree(root: TreeNode, split_node: TreeNode, outpath: Path):
    ts = TreeStyle()
    ts.mode = "c"
    # ts.show_leaf_name = False

    for n in root.traverse():
        if n == split_node:
            n.set_style(NodeStyle(fgcolor="red", size=100))
            subtree_colours = ["blue", "green", "purple"]
            for i, subtree in enumerate(n.children):
                for child_node in subtree.traverse():
                    child_node.set_style(NodeStyle(fgcolor=subtree_colours[i], size=20))
            break
        else:
            n.set_style(NodeStyle(fgcolor="black", size=0))

    root.render(str(outpath), w=180, units="mm", tree_style=ts)


if __name__ == "__main__":
    main()
