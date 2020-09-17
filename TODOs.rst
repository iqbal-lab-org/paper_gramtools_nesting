docs
======

[] section on data reuse. 
  --> Make available cortex pf3k vcfs, comas assemblies, clockwork (cortex) 1017 (1000 + 17 Comas) vcfs
  --> Accessions for 700 pf3k BAMs


reproducibility
================

[x] container: install python env from requirements.txt. See martin's clockwork container setup for how to interact with repo on local filesystem.
[x] container: RScript installing all dependencies in one go, called by singu def file.
[] Make sure container can run scripts which rely on nesting_paper py package (installed in dev mode from analysis/scripts)
[] Re-run all workflows + add Michael's singularity bind options

tb_bigdel
=========

[] analyse low-precision gramtools calls (esp. lower prec than graphtyper2)
    --> In evaluate_jvcfs.py script, instead of strictly looking at right end boundary for --region, extend up to the expected number of sites in the truth jvcf
    --> It seems I can get all nested calls correct but not the parent call correct :O; seen in N0004 tb_bigdel_28

msps_dimorphism
================

[] Trees: could Add deduplication step when catting all fa sequences together. RAxML complains about duplicates in the msa otherwise. It does produce a deduplicated file as well by default, disabled using `--no-seq-check`.



