analysis prerequisites
=========================

TB Comas data
--------------

This is ilmn reads and pacb assemblies of 17 samples. The public accessions are listed in analysis/input_data/mtuberculosis/Comas

Clockwork vcfs from ilmn reads against H37Rv reference are also listed in same directory. These are used to produce the initial set of deletion regions to construct prgs in.

NOTE: assemblies deemed improvable by martin: for eg, N1202 is missing the region approximately 899000 - 974500 from his analyses.

TB clockwork vcfs
------------------

Start with 13k ENA public run accessions from Mykrobe paper: https://figshare.com/articles/dataset/sample_data_tsv/7556789. 

These are searched against a database dump of clockwork-run var calling against all ENA TB samples: analysis/input_data/mtuberculosis/clockwork/clockwork_ena*.tsv.

Retrieved 6845 matches, giving links to the reads (clockwork_data_mykrobe_samples.tsv) and the vcfs (clockwork_vcfs_mykrobe_samples.tsv) on EBI's yoda cluster.

Note the vcfs.tsv points to directories on the cluster with subdirectories cortex, minos, samtools. Probably use cortex as avoids explaining/justifying minos using gramtools.

Ref genome
----------


results
========

vg
----

On vgteam/vg Gitter, 08/08/2020::

    Brice Letcher @bricoletc 15:53
    hello, is there a way in vg call to obtain calls corresponding to each snarl in the graph, including REF and null calls? i don't want to restrict to calls from an input vcf (-v) nor do I want to have only calls corresponding to ALTs, as I'm looking at a joint genotyping of several samples scenario

    Glenn Hickey @glennhickey 16:01
    @bricoletc No, you only get that with -v. This is a reasonable thing to want, and others have requested it. I will add an option asap.
