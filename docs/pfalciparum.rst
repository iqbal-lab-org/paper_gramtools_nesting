Describes some practical notes related to pfalciparum prgs and datasets.

analysis prerequisites
========================

input data
------------

prg construction
``````````````````

* pf3k vcfs: /nfs/research1/zi/mirror/oxford/cycloid/data3/projects/pf3k/discovery_gt/release3_for_first_paper/release3_for_first_paper on noah
* Ref genome and annotations: ftp://ftp.sanger.ac.uk/pub/project/pathogens/gff3/2015-08/


PACB/ILMN pairs
````````````````
* PACB assembled genomes: ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/
 --> What is Fort Lauderdale data usage policy?


* Corresp table between ILMN reads sample names and PACB assembled genome sample names. ENA accessions can be found in this paper: `Long read assemblies of geographically dispersed Plasmodium
falciparum isolates reveal highly structured subtelomeres`

* The data is downloaded via analysis/input_data/pfalciparum/pacb_ilmn/dl_*.sh scripts. Note I subsampled the reads for PfIT to 50% after dl using `seqtk` as it has > 2x more reads than others.

read sets
``````````

pf3k read sets:  /nfs/research1/zi/mirror/oxford/cycloid/data3/projects/pf3k/BAMs_release3.27.03.15

There are ~706 there. Could add more.


Questions on input data
------------------------

* Zam filtered out indels > 20bp to produce /nfs/research1/zi/bletcher/projects/Pf_benchmark/data/pf3k_and_DBMSPS1and2 ; why?
* pf3k read sets: how different are bams from same sample from one pf3k (or community project, or pf6) release to the next?



genes
======

* (R Amato suggested) EBA175, HRP2/3
* Miles et al 2016: 
  - DBLMSP 1 and DBLMSP2. These are part of the msp3 family.
  - AMA1
  - SERA5 
  - SURF1.2, 4.1, 8.2, 13.1, 14.1
  - MSP1, 3, 6
* Barry and Arnott 2014:

* Barry, Schultz, Buckee, Reeder 2009:
  [Pre-erythrocytic]
  - csp
  - trap
  - lsa1
  - glurp
  [Blood-stage]
  - eba175
  - ama1
  - msp1, 3, 4
  [Sexual-stage]
  - Pfs 48/45
* Malaria WHO Rainbow table
   [Pre-erythrocytic]
   - TRAP + ME epitopes (CS, LSA1, LSA3, STARP, EXP1, pb9)
   [Blood-stage]
   - MSP3 (with GLURP in some projects)
   - SERA5 N-terminal repeat domain
   - RH5
   - var2csa
   - AMA1-DiCo
   - P27A
   [Sexual stage]
   - Pfs25 

List of 4 genes to start with:
eba175, msp3.4, msp3.8, AMA1
