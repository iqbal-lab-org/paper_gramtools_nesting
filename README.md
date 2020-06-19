Workflows
===========

nocond_simulations
-------------------

This is to simulate paths in a given prg and genotype them, without conditions.
It can be used to evaluate changes in performance when changing, for eg, the genotyping model.

Plasmodium DBLMSPs
```````````````````

I re-use the prgs constructed in pf_surfants repository, on cluster.
I pulled results from `starting_prg` workflow in `pf_surfants`, which are per-gene prgs and non variant binary portions.
I used the 5kbp gene flank full.bed in `outputs/starting_prg/gene_list_1/beds` and cut it down to DBLMSP1 and 2 and the non var in between.
I call `pf_surfants/scripts/starting_prg/concat_prg.py` using that bed file to produce prgs. This is required so variant site numbers get rescaled across prgs.

I extract a `ref.fa` for that portion using samtools faidx on Pf ref with the region spanning the DBLMSPs noted in full.bed.

Evaluating genotyping
``````````````````````

Commits:

* 00eeef18 : no GCP
* 71e22c47 : first introduction of GCP. 


nestedness_simulations
-----------------------

This is to simulate paths through a nested (resp. non-nested prg), thread them through a non-nested (resp. nested) prg,
and compare performance across that condition change.

Plasmodium DBLMSPs
```````````````````

Prg production as per nocond_simulations.

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

Practical notes
````````````````

* Watch out in ref.fa, which is constructed from full.bed and Plasmo 3D7 ref in input_data/plasmodium/DBLMSPs, by default using this will name the CHROM in the vcf the gene name + sliced coordinates. Make sure there is a space between 'Pf3D7_10' and the sliced coords.

* input_data/plasmodium/prgs/DBLMSPs has a full.bed which describes the prg/ref.fa relative to 3D7 and a genes_in_prg.bed describing the genes only (not with the 5k flanks as in full.bed) relative to the prg. I checked that the sequences are the same.

TODOs
```````

[] Possible rmduplicates step when catting all fa sequences together. RAxML complains about duplicates in the msa otherwise. It does produce a deduplicated file as well by default, disabled using `--no-seq-check`.
