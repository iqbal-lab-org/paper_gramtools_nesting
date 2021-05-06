from typing import List

def get_samples(sample_tsv: str) -> List[str]:
    samples = []
    with open(sample_tsv) as fin:
        for line in fin:
            tabs = line.split("\t")
            if len(tabs) < 3:
                continue
            to_process = tabs[2]
            if to_process == "1":
                samples.append(tabs[0])
    return samples



def get_tree_genes(bed_path: str) -> List[str]:
    res = list()
    with open(bed_path) as fin:
        for line in fin:
            entry_name = line.strip().split("\t")[3]
            res.append(entry_name)
    return res
