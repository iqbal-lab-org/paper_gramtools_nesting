set -eu
base_dir="/nfs/leia/research/iqbal/bletcher/analyses/nesting_paper"
samp_accessions="${base_dir}/analysis/input_data/pfalciparum/pacb_ilmn/data_accessions.tsv"

output_dir="${base_dir}/analysis/input_data/pfalciparum/pacb_ilmn/ilmn_reads"
mkdir -p "$output_dir"

log_dir="${base_dir}/analysis/logs/downloads"
mkdir -p "$log_dir"

IFS=$'\n'
for line in $(cat ${samp_accessions})
do
	if [[ ${line:0:1} == '#' ]] || [[ ${line:0:6} == 'sample' ]] ; then continue; fi
	echo $line
	IFS=$'\t'; elems=($line)
	ilmn_read_accession="${elems[2]}"
	samp_name="${elems[0]}"
	bsub.py 1 \"${samp_name}_ena_dl\" -o ${log_dir}/${samp_name}_ilmn_reads.o -e ${log_dir}/${samp_name}_ilmn_reads.e "enaDataGet -f fastq -d ${output_dir}/${samp_name} $ilmn_read_accession"
done

