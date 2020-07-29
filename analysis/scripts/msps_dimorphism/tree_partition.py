import sys
from pathlib import Path
from typing import List, NamedTuple

from ete3 import Tree, TreeNode, TreeStyle, NodeStyle


def usage(msg: str = ""):
    if msg != "":
        print(f"ERROR: {msg}")
    print(
        f"{sys.argv[0]} newick_file output_file \n Outputs a basal partition of the tree, with each partition as a list of node names on a line"
    )
    exit(1)


NodeNames = List[str]


class PartitionResult(NamedTuple):
    groups: List[NodeNames]
    split_node: TreeNode


def get_partition(tree: Tree) -> PartitionResult:
    cur_node = tree
    groups = list()
    while True:
        no_leaf_children = [child for child in cur_node.children if not child.is_leaf()]
        if len(no_leaf_children) == 1:
            cur_node = no_leaf_children[0]
        elif len(no_leaf_children) == 0:
            raise ValueError("Could not find a partition")
        else:
            largest_size = -1
            largest_subtree = None
            total_size = 0
            for subtree in no_leaf_children:
                its_size = len(subtree)
                if its_size > largest_size:
                    largest_size = its_size
                    largest_subtree = subtree
                total_size += its_size
            if largest_size / total_size > 0.9:
                cur_node = largest_subtree
            else:
                for subtree in no_leaf_children:
                    its_nodes = [
                        node.name for node in subtree.traverse() if node.is_leaf()
                    ]
                    groups.append(its_nodes)
                break
    return PartitionResult(groups, cur_node)


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


def main():
    if len(sys.argv) != 3:
        usage()

    input_file = Path(sys.argv[1])
    if not input_file.exists():
        usage(f"{input_file} not found")

    output_prefix = Path(sys.argv[2]).resolve()
    if not output_prefix.parent.exists():
        output_prefix.parent.mkdir(parent=True)

    tree = Tree(str(input_file))
    partition = get_partition(tree)
    print(f"Found partition of sizes {[len(l) for l in partition.groups]}")

    with Path(f"{output_prefix}.basal_split").open("w") as fout:
        for l in partition.groups:
            fout.write("\t".join(l) + "\n")

    draw_tree(tree, partition.split_node, Path(f"{output_prefix}.pdf"))


if __name__ == "__main__":
    main()
