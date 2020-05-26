import pandas as pd

datasets = pd.read_csv(config["datasets"], sep="\t")
datasets = dict(zip(datasets["name"],datasets["base_path"]))

def get_data_path(wildcards):
    return datasets[wildcards.dataset]
