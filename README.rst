This is the repository for the gramtools nested variation paper analysis. It is divided into snakemake workflows, one folder each in analysis/workflows.


Workflows
===========

Configuration is always done via the workflow yaml file in analysis/configs. If a workflow can use one of several yaml (eg make_prgs), the desired one needs to be used in the workflow's Snakefile.

make_prgs
----------
Takes a file with a gene list + a .gff (TODO: OR a bed file of genomic regions) and builds a prg of those regions from a fasta ref + initial vcfs.

Can parameterise `make_prg` with 'max_nesting' and 'min_match_length' to make a set of different prgs for a given region list.


pacb_ilmn_validation
---------------------

Takes a PRG, truth assemblies from samples and ilmn reads from the same samples. Makes calls of the samples against the prg/ref genome and validates calls made in regions interest.

Truth assemblies and ilmn read sets should be downloaded before running the workflow. They should be placed in analysis/input_data/{dataset}/pacb_ilmn. For eg for pfalciparum I have download scripts in analysis/input_data/pfalciparum/pacb_ilmn/dl_*.sh

pfalciparum
````````````

I listed the ENA accessions for the ILMN and PACB reads in input/pfalciparum/pacb_ilmn/data_accessions.tsv. Note that when looking for the ILMN sample accessions in the pf3k release 3 metadata file, I find none. The release 3 data is what was used to get variants by cortex and thus to make the starting_prg. In pf3k release 5 metadata, I do find 3 accessions: ERS740936 (KH02), ERS740937(KE01) and ERS740940(GN01). 



nocond_simulations
-------------------

This is to simulate paths in a given prg and genotype them, without comparing conditions (such as PRG nestedness).
It can be used to evaluate changes in performance when changing, for eg, the genotyping model.


Evaluating genotyping
``````````````````````

Commits:

* 00eeef18 : no GCP
* 71e22c47 : first introduction of GCP. 


nestedness_simulations
-----------------------

This is to simulate paths through a nested (resp. non-nested prg), thread them through a non-nested (resp. nested) prg,
and compare performance between the two conditions.

Plasmodium DBLMSPs
```````````````````

Here is the output of concat_prg.py on the non_nested (mn1_mml7) data:
```
INFO:root:Processing: DBLMSP
INFO:root:Cumulative len prg: 37304
INFO:root:Cumulative num sites: 454

INFO:root:Processing: nonvar_12
INFO:root:Cumulative len prg: 44508
INFO:root:Cumulative num sites: 454

INFO:root:Processing: DBLMSP2
INFO:root:Cumulative len prg: 100720
INFO:root:Cumulative num sites: 872
```

And for the nested equivalent (mn5_mml7):

```
INFO:root:Processing: DBLMSP
INFO:root:Cumulative len prg: 23676
INFO:root:Cumulative num sites: 815

INFO:root:Processing: nonvar_12
INFO:root:Cumulative len prg: 30880
INFO:root:Cumulative num sites: 815

INFO:root:Processing: DBLMSP2
INFO:root:Cumulative len prg: 68680
INFO:root:Cumulative num sites: 1717
```

--> Nested prg has 1.96 x more variant sites, and 0.68 x the num of characters!


msps_dimorphism
-----------------

This is to analyse dimorphisms in DBLMSP1 and DBLMSP2 from pf3k genotyped samples on the DBLMSP prg.



Running on cluster
====================

Requirements for running
``````````````````````````

* Snakemake==v5.14.0
* Singularity>=v3.4.0-1

Steps for running
```````````````````
* Requires singu container image in container/built. Can be built for example running `sudo singularity build container/built/singu.sif container/singu_def.def`. 

* Current working directory when running any workflow should be the git top-level directory(where this file is)

* On ebi cluster: `module load singularity/3.[45].[0-9]`

* Requires  lsf profile for snakemake: https://github.com/Snakemake-Profiles/snakemake-lsf. Requires as in cluster config yaml file is being deprecated.
The non-defaults configs are: LSF_UNIT_FOR_LIMITS=MB, default_cluster_logdir=run/logs/lsf_profile

* `bash analysis/cluster_submit.sh workflow_name`!

