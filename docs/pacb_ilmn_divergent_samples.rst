Exploring samples with high edit distance to pacb truths
=========================================================

First, these can be obtained like so::

    awk '{if (DBLMSP" && $3 == "gramtools_genotype"){print $0}}' stats.tsv  | sort -n -k4'

for gramtools_genotype on DBLMSP.

This shows four samples with edit distance > 4%: KH01, GN01, GA01, TG01.

Note that for **all four**, the edit distance between cortex run on the sample and 3D7 is **larger**. If the assembly ran fine, this suggests we are missing variation in the PRG.


Do ilmn reads agree with pacb assembly?
````````````````````````````````````````
[Work from tmp_work/ilmn_pacb_alignments]

I checked whether assembly and ilmn reads agreed by mapping latter onto former for PfGA01 (8.2% edit distance). Then take annotation of the assembly and find DBLMSP coordinates on it. Then do samtools stats on that region, and find pb mismatch rate of 8x10e-3 and almost no indels. bcftools mpileup + call finds no variants. ART shows good pileup in whole region. So assembly looks fine.

The differences concentrate in a small region near the beginning of the gene (can see this using ```cat PfGA01.bam | sam2pairwise```. 

Are there input samples, that went into the prg, that are close to the pacb assembly?
`````````````````````````````````````````````````````````````````````````````````````
Then I check if any of the input DBLMSP gene sequences to make_prg map closely to the assembly::

    fout=prg_seqs/DBLMSP_geneOnly.fa; >$fout; for id in  $(grep "^>" prg_seqs/DBLMSP.fa); do reg=${id/'>'/}:5000-7094; samtools faidx prg_seqs/DBLMSP.fa $reg >> $fout; done
    bowtie2 -x bowtie_indexes/PfGA01 -U prg_seqs/DBLMSP_geneOnly.fa -f --no-unal > matches.sam
    grep -E "NM:i:[[:digit:]]+" -o matches.sam | grep -E "[[:digit:]]+" -o | sort -n -r | uniq -c


And surprisingly, yes. PfGA01 has (2094 * 0.0821) = 171 differences while there are 7 samples with only 70 differences in the set of sequences that go into msa and then make_prg. 

So msa or make_prg construction means we lose some variants, or gramtools cannot find these variants in the graph

Additionally, note cortex does not find any more variants when run against the personalised reference that it 0.0821 diverged. 


Is cortex the issue and does samtools do better?
````````````````````````````````````````````````
So I run samtools var calling (in tmp_work/ilmn_pacb_alignments/run_samtools). I get the gramtools personalised ref and do the following::

    samtools view -h PfGA01.bam PfGA01_10:1428143-1430242 > run_samtools/DBLMSP_PfGA01.sam
    samtools fastq DBLMSP_PfGA01.sam > DBLMSP_PfGA01.fq
    bwa mem PfGBA01_gramtools_genotype_pers_ref.fa DBLMSP_PfGA01.fq > mapped.sam
    samtools sort -O BAM -o mapped.bam
    bcftools mpileup -f PfGBA01_gramtools_genotype_pers_ref.fa mapped.bam | bcftools call -mv -Ob -o calls.bcf
    bcftools index calls.bcf
    samtools faidx PfGBA01_gramtools_genotype_pers_ref.fa Pf3D7_10_v3:1413200-1415294 | bcftools consensus calls.vcf -H 1 > DBLMSP_withvars.fa

The first line uses the coordinates of DBLMSP on the PfGA01 pacb assembly.

This reduced the edit distance by 10 (when mapping the .fa to the pacb assembly with bowtie2), so makes almost no difference.
    

Back to low divergence input samples: can we obtain them in the graph?
````````````````````````````````````````````````````````````````````````````

Take one sample with edit distance 70 to the pacb assembly: PH0942-Cx. Simulate reads from it with art_ilmn, and genotype them on DBLMSP mn5_mml7 prg (gram build in nestedness workflow outputs).

The resulting personalised ref has edit distance 1 to this sample. So it can be found (bar this annoying 1).


