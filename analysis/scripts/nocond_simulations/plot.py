import sys
from pathlib import Path

import pandas as pd

def usage():
    print(f"Usage: {sys.argv[0]} input_tsv output_dir")
    exit(1)

def get_correctness(dataframe_row) -> str:
    series = dataframe_row[1]
    if series.res_has_call:
        if series.res_is_correct:
            return "TP"
        else:
            return "FP"
    else:
        if series.res_is_correct:
            return "TN"
        else:
            return "FN"

def main():
    if len(sys.argv) != 3:
        usage()

    input_tsv = Path(sys.argv[1]).resolve()
    if not input_tsv.exists():
        print(f"{input_tsv} not found")
        usage()

    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(input_tsv, dep="\t")
    correctness = list(map(get_correctness, data.iterrows()))
    data["correctness"] = correctness

    
