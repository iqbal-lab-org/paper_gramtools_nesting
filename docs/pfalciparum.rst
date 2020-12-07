Describes some practical notes related to pfalciparum prgs and datasets.

analysis prerequisites
========================

input data
------------

prg construction
``````````````````

* pf3k vcfs: /nfs/research1/zi/mirror/oxford/cycloid/data3/projects/pf3k/discovery_gt/release3_for_first_paper/release3_for_first_paper on noah
  --> Distributed on Zenodo
* Ref genome and annotations: ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/2015-08/


PACB/ILMN pairs
````````````````
* PACB assembled genomes: ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/
  Annotation files are also available, under Isolates/.


* Corresp table between ILMN reads sample names and PACB assembled genome sample names. ENA accessions can be found in this paper: `Long read assemblies of geographically dispersed Plasmodium
falciparum isolates reveal highly structured subtelomeres`

* The data is downloaded via analysis/input_data/download_data/pf_validation_dl_*.sh scripts. Note I subsampled the reads for PfIT to 50% after dl using `seqtk` as it has > 2x more reads than others.

read sets
``````````

pf3k read sets:  /nfs/research1/zi/mirror/oxford/cycloid/data3/projects/pf3k/BAMs_release3.27.03.15
There are 706 there. Could add more.
    --> ENA read accessions used for reproducibility



Questions on input data
------------------------

* Zam filtered out indels > 20bp to produce /nfs/research1/zi/bletcher/projects/Pf_benchmark/data/pf3k_and_DBMSPS1and2
* pf3k read sets: how different are bams from same sample from one pf3k (or community project, or pf6) release to the next?

