from typing import List

def get_bams(filelist_path: str) -> List[str]:
        res = list()
	with open(filelist_path) as fin:
            for line in fin:
                assert(line.rstrip().endswith(".bam"))
                res.append(line.rstrip().replace(".bam",""))
        return res
