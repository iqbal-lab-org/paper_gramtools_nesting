Intro
======
This is the repository for the paper *Enabling multiscale variation analysis with genome graphs*. 

It is divided into snakemake workflows, one folder each in analysis/workflows.

Setup
=======

Requirements for setup
--------------------------

* Python >=3.6
* [Singularity>=v3.5.0](https://sylabs.io/guides/3.5/user-guide/>) (Recommended: 3.7.0)


Input data
------------

Here we list all the datasets used; for each, either public accessions are stored in a file or the data is available for download on Zenodo. 
We provide scripts and commands to obtain all the data, see [below](#Steps-for-setup).

Data                      | Access
----                      | -----
*M. tuberculosis* vcfs of 1,017 samples | Zenodo
*P. falciparum* vcfs of 2,498 samples | Zenodo
*M. tuberculosis* validation PacBio assemblies of 17 samples | Zenodo
*M. tuberculosis* validation Illumina read sets of 17 samples | ENA accessions listed in `analysis/input_data/mtuberculosis/pacb_ilmn/data_accessions.tsv`
*P. falciparum* 14 Illumina read sets and matched PacBio assemblies | ENA accessions listed in `analysis/input_data/pfalciparum/pacb_ilmn/data_accessions.tsv`
*P. falciparum* 706 Illumina read sets from Pf3k release 5 | ENA accessions listed in `analysis/input_data/pfalciparum/pf3k/pf3k_release_5.tsv`


Steps for setup
-------------------

All commands work only when issued from the root of this project.

Setup a python virtual environment and obtain snakemake:
    
```sh
    git clone https://github.com/iqbal-lab-org/paper_gramtools_nesting && cd paper_gramtools_nesting
    python3 -m venv venv && . venv/bin/activate
    pip3 install pip==20.0.2 
    pip3 install -r pyrequirements.txt
```

Obtain the singularity container:

```sh
    # Download container:
    mkdir -p container/built && wget https://zenodo.org/record/4147302/files/1_container.sif?download=1 -O container/built/singu.sif
    # It is possible to build the container directly; this should take ~30 minutes:
    # sudo singularity build --no-cleanup container/built/singu.sif container/singu_def.def 
```

Obtain the input data:

```sh
    # Some data download uses enaBrowserTools and seqtk, which are installed in container.
    singu_command="singularity exec container/built/singu.sif"

    ## P. falciparum data ##
    "$singu_command" bash analysis/input_data/download_data/pf_dl_ilmn_ena.sh
    "$singu_command" bash analysis/input_data/download_data/pf_dl_pacb_assemblies.sh

    # Below downloads 706 fastqs of ~1GB each; I recommend modifying the script to submit in parallel to a cluster (adding in singularity command too). Downloads from ENA server, several reruns may be required if server throws any error.
    bash analysis/input_data/download_data/pf3k_dl_ilmn_all.sh

    pf_vcfs="analysis/input_data/pfalciparum/pf3k/vcfs"
    mkdir -p "$pf_vcfs"
    wget https://zenodo.org/record/4147302/files/3_pf_vcfs.tar?download=1  -O 3_pf_vcfs.tar
    tar -xf 3_pf_vcfs.tar

    ## M. tuberculosis data ##
    tb_vcfs="analysis/input_data/mtuberculosis/clockwork/vcfs"
    mkdir -p "$tb_vcfs"
    
    wget https://zenodo.org/record/4147302/files/2_mtb_vcfs.tar?download=1 -O 2_mtb_vcfs.tar
    tar -xf 2_mtb_vcfs.tar
    wget https://zenodo.org/record/4147302/files/4_mtb_validation_assemblies.tar?download=1 -O 4_mtb_assemblies.tar
    tar -xf 4_mtb_assemblies.tar

    "$singu_command" python3 analysis/input_data/download_data/mtb_dl_ilmn_ena.py
```


How to run a workflow
----------------------

```sh
    . venv/bin/activate
    WORKFLOW=<workflow_name>
    # Add -nq to the command below to only preview what would be run.
    snakemake -s analysis/workflows/"${WORKFLOW}"/Snakefile --use-singularity --verbose
```

We have run the analyses on an LSF Cluster, executing the following:

* Install the snakemake lsf profile https://github.com/Snakemake-Profiles/snakemake-lsf;
  the non-defaults configs are: LSF_UNIT_FOR_LIMITS=MB, default_cluster_logdir=analysis/logs/lsf_profile
* Run submission script:

```sh
    . venv/bin/activate
    bash analysis/cluster_submit.sh <workflow_name>
```

The submission script contains specific filepaths on our cluster which need modifying.


Workflow details
============================

make_prgs
-----------

No prerequisites.

In analysis/workflows/make_prgs, at top of Snakefile, there are two `configfile:` directives, pointing to `pfalciparum.yaml` and `mtuberculosis.yaml`. Uncomment the line corresponding to which analysis to run.

Now run :

```sh
    bash analysis/cluster_submit.sh make_prgs
```

This produces inputs required by `nestedness_simulations`, `pacb_ilmn_validation`, `msps_dimorphism` and `tb_bigdel` workflows.

nestedness_simulations
-----------------------

Requires: [make_prgs](#make_prgs) on pfalciparum.yaml

This produces the results of the section `Validation of nested genotyping with simulated data`. Because paths in the graphs, and reads from the paths, are randomly simulated, the exact same result will not be produced. I ran the workflow multiple times and got the same performance each time.

pacb_ilmn_validation
---------------------
Requires: [make_prgs](#make_prgs) on pfalciparum.yaml

This produces the results of the section `Validation of nested genotyping with simulated data`.

pacb_ilmn_prg_closest_pf
-------------------------
Requires: pacb_ilmn_validation

This produces the results of the section: `Comparing gramtools personalised reference with best of reference panel`.


tb_bigdel
----------
Requires: make_prgs on mtuberculosis.yaml ; vg_make_prgs

This produces the results of the section: `Application: unified SNP and large deletion analysis in M. tuberculosis`.

msps_dimorphism
-----------------
Requires: make_prgs on pfalciparum.yaml

This produces the results of the section: `Application: charting SNPs on top of alternate haplotypes`.


Development
============

While working on the paper, it is not convenient as you have to rebuild and transfer the whole container image if anything changes (eg one line in gramtools). For development, it is best to work with a python `venv` at top-level. Then `pip install -r pyrequirements.txt` and `pip install -r container/pyrequirements.txt`. Then locate gramtools on cluster, cmake/make it, and `pip install -e` it from inside the `venv`. And in the file `analysis/cluster_submit.sh`, remove the `--use-singularity` argument. All tools must be available on cluster.

