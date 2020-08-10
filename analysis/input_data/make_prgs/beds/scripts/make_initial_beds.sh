#pfalciparum: make bed from gene_list + gff of genome
bash gff_to_bed.sh ../gene_lists/pf_DBLMSPs.txt ../../pfalciparum/ref_genome/Pfalciparum.gff3 pf_DBLMSPs.bed
bash gff_to_bed.sh ../gene_lists/pf_4surfants.txt ../../pfalciparum/ref_genome/Pfalciparum.gff3 pf_4surfants.bed

# mtuberculosis: make bed from file of vcf file names and deletion size
bash bigdel_vcfs_to_bed.sh ../../mtuberculosis/clockwork/vcfs_Comas_pacb_ilmn.txt 100 tb_bigdel.bed
