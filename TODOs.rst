docs
======

[] section on data reuse. 
  --> Make available cortex pf3k vcfs, comas assemblies, clockwork (cortex) 1017 (1000 + 17 Comas) vcfs
  --> Accessions for 700 pf3k BAMs


reproducibility
================

[x] container: install python env from requirements.txt. See martin's clockwork container setup for how to interact with repo on local filesystem.
[x] container: RScript installing all dependencies in one go, called by singu def file.
[x] container: add enaDataGet for downloading from ENA?
[x] container: ensure can run scripts which rely on nesting_paper py package (installed in dev mode from analysis/scripts)
[] add Michael's singularity bind options + check singularity workdir

[x] Re-run all workflows 
[] pf3k bams: script to download them

tb_bigdel
=========

[x] analyse low-precision gramtools calls (esp. lower prec than graphtyper2)
    [x] In evaluate_jvcfs.py script, instead of strictly looking at right end boundary for --region, extend up to the expected number of sites in the truth jvcf
    [x] I could get all nested calls correct but not the parent call correct, fixed by improving gramtools gtyping model 

msps_dimorphism
================

[] Trees: could Add deduplication step when catting all fa sequences together. RAxML complains about duplicates in the msa otherwise. It does produce a deduplicated file as well by default, disabled using `--no-seq-check`.



