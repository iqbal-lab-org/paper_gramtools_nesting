#!/usr/bin/env bash
args=("$@")
fasta_ref="${args[0]}"
reads_files="${args[1]}"
output_vcf="${args[2]}"
threads="${args[3]}"
cp $fasta_ref ref
bwa index ref
bwa mem -t $threads ref $reads_files -o mapped.sam
samtools sort -o mapped.bam -O BAM mapped.sam
samtools index mapped.bam
bcftools mpileup -f ref mapped.bam | bcftools call -O z -o $output_vcf -vm --ploidy 1
bcftools index $output_vcf
