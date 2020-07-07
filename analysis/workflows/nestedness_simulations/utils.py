import pandas as pd


input_table = pd.read_csv(config["datasets"], sep="\t")
datasets = set(input_table["name"])
conditions = [
    "nonested",
    "nested",
]  # Ordered, so that simulate paths from non-nested prg and induce them in nested prg
assert set(conditions) == set(input_table["nesting"])
dataset_access = dict()
for _, row in input_table.iterrows():
    dataset_access[f'{row["name"]}_{row["nesting"]}'] = {
        "prg": row["base_path"],
        "genome_ref": row["genome_ref"],
        "coords": row["coords"],
    }


def get_prg_path(wildcards):
    return dataset_access[f"{wildcards.dataset}_{wildcards.nesting}"]["prg"]


def get_genome_path(wildcards):
    return dataset_access[f"{wildcards.dataset}_{wildcards.nesting}"]["genome_ref"]


def get_coords(wildcards):
    return dataset_access[f"{wildcards.dataset}_{wildcards.nesting}"]["coords"]
