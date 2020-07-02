
make_prgs
=========

* Support for a bed input. Currently I use a list of genes + a gff3 to produce the bed in rule 'make_beds'


msps_dimorphism
================

* Draw trees with the basal split node coloured, to relate the groups to the networks showing divergence
Thus probably produce the basal splits in plot_tree.R directly, deprecating the python script that does it.


* Add the required r packages in the scripts that need them, perhaps as a call to biocManager::install / install.packages.
Or make an RScript installing all dependencies in one go, called by singu def file.


nestedness_simulations
========================

* downloaded art_illumina binary fails looking for 'libgsl.so.0' on singu image. Maybe compile it from source.
