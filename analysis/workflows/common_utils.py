# Get gramtools commit
gmtools_commit_script = Path(config["scripts"]) / "gmtools_commit.py"
GMTOOLS_COMMIT = sp_run(
    ["singularity", "exec", config["container"], "python3", gmtools_commit_script],
    stdout=PIPE,
    universal_newlines=True,
)
GMTOOLS_COMMIT.check_returncode()
GMTOOLS_COMMIT = GMTOOLS_COMMIT.stdout.strip()
