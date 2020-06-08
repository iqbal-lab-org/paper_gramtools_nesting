import pandas as pd

from gramtools import version

datasets = pd.read_csv(config["datasets"], sep="\t")
datasets = dict(zip(datasets["name"],datasets["base_path"]))

def get_data_path(wildcards):
    return datasets[wildcards.dataset]

_ , version_dict = version.report()
GMTOOLS_COMMIT = version_dict.get("last_git_commit_hash")
if GMTOOLS_COMMIT is None:
    raise ValueError("Could not get gramtools commit hash from gramtools version. Gramtools is probably not compiled.")
else:
    GMTOOLS_COMMIT = GMTOOLS_COMMIT[:8]
