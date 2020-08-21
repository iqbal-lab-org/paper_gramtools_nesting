# This script maps ilmn reads to i)ref genome ii)assembly for samples for which column 4, which indicates whether a deletion has been recovered, is set to 0.
# It also maps the reference portion of the genome around the deletion to the assembly to obtain its coordinates in the assembly.
# The outputs allow manual inspection of said deletions using ACT (Artemis Comparison Tool)

reads_dir="analysis/input_data/mtuberculosis/pacb_ilmn/ilmn_reads"
assembly_dir="analysis/outputs/tb_bigdel/bowtie_indexes"
genome="analysis/input_data/mtuberculosis/ref_genome/H37Rv.fasta"
genome_name="NC_000962.3"
minimap_eval=""

base_dir="analysis/outputs/tb_bigdel/alignments/"
ref_portions="${base_dir}/ref_portions" && mkdir -p $ref_portions
ref_idx="${base_dir}/ref_idx" && mkdir -p $ref_idx
aln_asm="${base_dir}/ilmn_alignments_asm" && mkdir -p $aln_asm
aln_ref="${base_dir}/ilmn_alignments_ref" && mkdir -p $aln_ref

flank_size=3000
echo "flank size: $flank_size"

if [[ ! -e ${ref_idx}/H37Rv.1.bt2 ]]; then
	bowtie2-build $genome ${ref_idx}/H37Rv
fi

IFS=$'\n'
for line in $(awk '{if ($4 == 0) print $0}' $minimap_eval); 
do
	IFS=$'\t';fields=($line)
	start=$((${fields[0]} + 1))
	pre_start=$((start - flank_size))
	stop=$((start + ${fields[1]} + flank_size))
	sample=${fields[2]}
	fname=${sample}_${start}_${fields[1]}
	samtools faidx $genome ${genome_name}:${pre_start}-${stop} > ${ref_portions}/${fname}.fa
	res="ref_portions/${fname}.sam"
	if [[ ! -e $res ]];then
	bowtie2 -x ${assembly_dir}/${sample} -U ${ref_portions}/${fname}.fa -f > $res
	fi

	reads_1="${reads_dir}/${sample}/reads_1.fastq.gz"
	reads_2="${reads_dir}/${sample}/reads_2.fastq.gz"

	res_sam="${ilmn_alignments_asm}/${fname}.sam"
	res_bam="${ilmn_alignments_asm}/${fname}.bam"
	if [[ ! -e $res_bam ]];then
		bsub.py --threads 4 2 "${fname}" "\"bowtie2 -x ${assembly_dir}/${sample} -1 ${reads_1} -2 ${reads_2} -S $res_sam -p 4 && samtools sort $res_sam -O BAM -o $res_bam\""
	fi

	res_sam="${ilmn_alignments_ref}/${fname}.sam"
	res_bam="${ilmn_alignments_ref}/${fname}.bam"
	if [[ ! -e $res_bam ]];then
		bsub.py --threads 4 2 "${fname}" "\"bowtie2 -x ${ref_idx}/H37Rv -1 ${reads_1} -2 ${reads_2} -S $res_sam -p 4 && samtools sort $res_sam -O BAM -o $res_bam\""
	fi

done
