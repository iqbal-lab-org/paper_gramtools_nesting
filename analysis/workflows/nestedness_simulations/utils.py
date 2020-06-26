import pandas as pd


input_table = pd.read_csv(config["datasets"], sep="\t")
datasets = set(input_table["name"])
conditions = [
    "nonested",
    "nested",
]  # Ordered, so that simulate paths from non-nested prg and induce them in nested prg
assert set(conditions) == set(input_table["nesting"])
path_access = dict(
    zip(input_table["name"] + "_" + input_table["nesting"], input_table["base_path"])
)


def get_data_path(wildcards):
    return path_access[wildcards.dataset + "_" + wildcards.nesting]
