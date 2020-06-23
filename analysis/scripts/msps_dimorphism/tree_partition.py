import sys
from pathlib import Path
from typing import List

from ete3 import Tree


def usage(msg: str = ""):
    if msg != "":
        print(f"ERROR: {msg}")
    print(
        f"{sys.argv[0]} newick_file output_file \n Outputs a basal partition of the tree, with each partition as a list of node names on a line"
    )
    exit(1)


Nodes = List[str]


def get_partition(tree: Tree) -> List[Nodes]:
    cur_node = tree
    result = list()
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
                    result.append(its_nodes)
                break
    return result


def main():
    if len(sys.argv) != 3:
        usage()

    input_file = Path(sys.argv[1])
    if not input_file.exists():
        usage(f"{input_file} not found")

    output_file = Path(sys.argv[2]).resolve()
    if not output_file.parent.exists():
        output_file.parent.mkdir(parent=True)

    tree = Tree(str(input_file))
    partition = get_partition(tree)
    print(f"Found partition of sizes {[len(l) for l in partition]}")

    with output_file.open("w") as fout:
        for l in partition:
            fout.write("\t".join(l) + "\n")


if __name__ == "__main__":
    main()
