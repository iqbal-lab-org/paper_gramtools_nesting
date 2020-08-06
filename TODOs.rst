docs
======

* section on data reuse. 
  -->TB vcfs and bams, which samples can be used freely? and Comas samples?


reproducibility
================

* container: install python env from requirements.txt. See martin's clockwork container setup for how to interact with repo on local filesystem.
* R scripts: install all required r packages perhaps as a call to biocManager::install / install.packages, in each script.
Or make an RScript installing all dependencies in one go, called by singu def file.

tb_bigdel
=========

* make heatmap of Comas samples, in regions with bigdels only

msps_dimorphism
================

* Trees: could Add deduplication step when catting all fa sequences together. RAxML complains about duplicates in the msa otherwise. It does produce a deduplicated file as well by default, disabled using `--no-seq-check`.


nestedness_simulations
========================

