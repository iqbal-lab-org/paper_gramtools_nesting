This is the repository for the gramtools nested variation paper analysis. It is divided into snakemake workflows, one folder each in analysis/workflows.

They are designed to be run in a cluster environment with a singularity container containing all dependencies. 

Setup
=======

Requirements for setup
--------------------------

* Python >=3.6
* Singularity>=v3.4.0-1 + root (or `fakeroot <https://sylabs.io/guides/3.5/user-guide/fakeroot.html>`_) privilege for building container

Steps for setup
-------------------

Setup commands and run commands should always be issued from the root of this project.

Setup a python virtual environment and obtain snakemake::
    
    python3 -m venv venv && . venv/bin/activate
    pip3 install pip==20.0.2 
    pip3 install -r pyrequirements.txt

Build singularity container::

    sudo singularity build container/built/singu.sif container/singu_def.def 

Obtain input data::

    TODO

How to run a worfklow
----------------------
::

    . venv/bin/activate
    TODO: Below is ebi cluster + lsf-specific
    module load singularity/3.5.0
    bash analysis/cluster_submit.sh <workflow_name>

* On ebi cluster: `module load singularity/3.[45].[0-9]`

* Requires lsf profile for snakemake: https://github.com/Snakemake-Profiles/snakemake-lsf. 
    The non-defaults configs are: LSF_UNIT_FOR_LIMITS=MB, default_cluster_logdir=run/logs/lsf_profile


Order to run workflows in
============================

make_prgs
`````````
Requires: None
Run ::

   bash analysis/cluster_submit.sh make_prgs


In analysis/workflows/make_prgs, at top of Snakefile, there are two `configfile:` directives, pointing to `pfalciparum.yaml` and `mtuberculosis.yaml`. You just ran pfalciparum. Comment out that line, and uncomment the mtuberculosis one.  Run again ::
   bash analysis/cluster_submit.sh make_prgs

This produces inputs required by `nestedness_simulations`, `pacb_ilmn_validation`, `msps_dimorphism` and `tb_bigdel` workflows.

nestedness_simulations
```````````````````````
Requires: make_prgs on pfalciparum.yaml

This produces the simulation results of the paper. Because paths in the graphs, and reads from the paths, are randomly simulated, the exact same result will not be produced. I have rerun the workflow once more and confirmed the results in the paper.

pacb_ilmn_validation
`````````````````````
Requires: make_prgs on pfalciparum.yaml

This produces part of the genotyping in DBLMSP2 results of the paper: performance of gramtools compared to ref-based variant callers.


msps_dimorphism
`````````````````
Requires: make_prgs on pfalciparum.yaml

This produces the dimorphism in DBLMSP2 results of the paper.


tb_bigdel
``````````
Requires: make_prgs on mtuberculosis.yaml ; vg_make_prgs

This produces the benchmark results of the paper: analysis of large deletions and small variants under them in M. tuberculosis with comparison to vg and graphtyper2.


Development
============

While working on the paper, it is not convenient as you have to rebuild and transfer the whole container image if anything changes (eg one line in gramtools). For development, it is best to work with a python `venv` at top-level. Then `pip install -r pyrequirements.txt` and `pip install -r container/pyrequirements.txt`. Then locate gramtools on cluster, cmake/make it, and `pip install -e` it from inside the `venv`. And in the workflows, comment out `container:` line. All tools must be available on cluster.

