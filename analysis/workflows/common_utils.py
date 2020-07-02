from pathlib import Path
from subprocess import run as sp_run, PIPE

# Get gramtools commit through singularity container
gmtools_commit_script = Path(config["scripts"]) / "gmtools_commit.py"
GMTOOLS_COMMIT = sp_run(
    ["singularity", "exec", config["container"], "python3", gmtools_commit_script],
    stdout=PIPE,
    universal_newlines=True,
)
GMTOOLS_COMMIT.check_returncode()
GMTOOLS_COMMIT = GMTOOLS_COMMIT.stdout.strip()


def mk_output_dirs():
    """For each variable starting with 'output', makes the directory name it holds"""
    for variable in filter(lambda name: name.startswith("output"), dir()):
        Path(eval(variable)).mkdir(exist_ok=True, parents=True)
