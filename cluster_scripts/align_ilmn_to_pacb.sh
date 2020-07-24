module load singularity/3.5.0
output_dir="/hps/nobackup/research/zi/bletcher/nesting_paper/analysis/outputs/pacb_ilmn_validation/align_ilmn_to_pacb"
idx_dir="/hps/nobackup/research/zi/bletcher/nesting_paper/analysis/outputs/pacb_ilmn_validation/bowtie_indexes"
log_dir="${output_dir}/logs"
mkdir -p $output_dir $log_dir
for sample in PfGN01 PfKH01 PfGA01 PfTG01; do
	log=${log_dir}/${sample}
	reads=($(find /hps/nobackup/research/zi/bletcher/nesting_paper/analysis/input_data/pfalciparum/pacb_ilmn/ilmn_reads/${sample} -name "*.fastq.gz"))
	bsub.py 10 $sample -o ${log}.o -e ${log}.e "singularity exec /hps/nobackup/research/zi/bletcher/nesting_paper/container/built/singu.sif bowtie2 -x ${idx_dir}/${sample} -1 ${reads[0]} -2 ${reads[1]} -S ${output_dir}/${sample}.sam && samtools sort -O BAM -o ${output_dir}/${sample}.bam ${output_dir}/${sample}.sam && rm ${output_dir}/${sample}.sam"
done
