from pathlib import Path
from subprocess import run as sp_run, PIPE


def get_gmtools_commit():
    """Get gramtools commit version through singularity container it is installed in"""
    gmtools_commit_script = Path(config["scripts"]) / "gmtools_commit.py"
    GMTOOLS_COMMIT = sp_run(
        ["singularity", "exec", config["container"], "python3", gmtools_commit_script],
        stdout=PIPE,
        universal_newlines=True,
    )
    GMTOOLS_COMMIT.check_returncode()
    GMTOOLS_COMMIT = GMTOOLS_COMMIT.stdout.strip()
    if GMTOOLS_COMMIT is None:
        raise ValueError("There is no gramtools commit")
    return GMTOOLS_COMMIT


def mk_output_dirs(variables):
    """For each variable starting with 'output', makes the directory name it holds"""
    for variable in filter(lambda name: name.startswith("output"), variables):
        Path(eval(variable)).mkdir(exist_ok=True, parents=True)
