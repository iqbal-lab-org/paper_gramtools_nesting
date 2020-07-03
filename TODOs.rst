make_prgs
=========

* Support for a bed input. Currently I use a list of genes + a gff3 to produce the bed in rule 'make_beds'


pacb_ilmn_validation
======================

* Plot histograms rather than boxplots of edit distances

msps_dimorphism
================

* Add the required r packages in the scripts that need them, perhaps as a call to biocManager::install / install.packages.
Or make an RScript installing all dependencies in one go, called by singu def file.


* Draw trees with the basal split node coloured, to relate the groups to the networks showing divergence
Thus probably produce the basal splits in plot_tree.R directly, deprecating the python script that does it.


* Add deduplication step when catting all fa sequences together. RAxML complains about duplicates in the msa otherwise. It does produce a deduplicated file as well by default, disabled using `--no-seq-check`.


nestedness_simulations
========================

