docs
======

[] section on data reuse. 
  --> Make available cortex pf3k vcfs, comas assemblies, clockwork (cortex) 1017 (1000 + 17 Comas) vcfs
  --> Accessions for 700 pf3k BAMs


reproducibility
================

[x] container: install python env from requirements.txt. See martin's clockwork container setup for how to interact with repo on local filesystem.
[x] container: RScript installing all dependencies in one go, called by singu def file.
[] Re-run all workflows + add Michael's singularity bind options

tb_bigdel
=========

[] analyse low-precision gramtools calls (esp. lower prec than graphtyper2)

msps_dimorphism
================

[] Trees: could Add deduplication step when catting all fa sequences together. RAxML complains about duplicates in the msa otherwise. It does produce a deduplicated file as well by default, disabled using `--no-seq-check`.



