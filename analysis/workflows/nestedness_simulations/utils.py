import pandas as pd

from gramtools import version

input_table = pd.read_csv(config["datasets"], sep="\t")
datasets = set(input_table["name"])
conditions = list(set(input_table["nesting"]))
assert(len(conditions) == 2) 
path_access = dict(zip(input_table["name"] + "_" + input_table["nesting"],input_table["base_path"]))

def get_data_path(wildcards):
    return path_access[wildcards.dataset + "_" + wildcards.nesting]

_ , version_dict = version.report()
GMTOOLS_COMMIT = version_dict.get("last_git_commit_hash")
if GMTOOLS_COMMIT is None:
    raise ValueError("Could not get gramtools commit hash from gramtools version. Gramtools is probably not compiled.")
else:
    GMTOOLS_COMMIT = GMTOOLS_COMMIT[:8]
