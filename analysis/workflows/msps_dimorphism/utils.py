from typing import List

def get_bams(filelist_path: str) -> List[str]:
    res = list()
    with open(filelist_path) as fin:
        for line in fin:
            assert(line.rstrip().endswith(".bam"))
            res.append(line.rstrip().replace(".bam",""))
    return res

def get_tree_genes(bed_path: str) -> List[str]:
    res = list()
    with open(bed_path) as fin:
        for line in fin:
            entry_name = line.strip().split("\t")[3]
            res.append(entry_name)
    return res
