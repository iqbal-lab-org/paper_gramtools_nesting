from pathlib import Path
from subprocess import run as sp_run, PIPE


def get_gmtools_commit():
    """Get gramtools commit version through venv/singularity container it is installed in"""
    gmtools_commit_script = Path(config["scripts"]) / "gmtools_commit.py"
    gmtools_commit_script = str(gmtools_commit_script.resolve())
    try:  # Try in virtual environment first
        GMTOOLS_COMMIT = sp_run(
            ["python3", gmtools_commit_script],
            capture_output=True,
            text=True,
            check=True,
        )
    except:
        GMTOOLS_COMMIT = sp_run(
            [
                "singularity",
                "exec",
                str(Path(config["container"]).resolve()),
                "python3",
                gmtools_commit_script,
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    GMTOOLS_COMMIT = GMTOOLS_COMMIT.stdout.strip()
    return GMTOOLS_COMMIT


def mk_output_dirs(variables):
    """For each variable starting with 'output', makes the directory name it holds"""
    for variable in filter(lambda name: name.startswith("output"), variables):
        Path(eval(variable)).mkdir(exist_ok=True, parents=True)
