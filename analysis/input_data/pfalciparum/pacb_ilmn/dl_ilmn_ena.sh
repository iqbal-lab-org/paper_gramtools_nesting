set -eu

[[ $# != 1 ]] && echo "usage: $0 base_directory (root path of the repository)" && exit 1
base_dir=$(realpath $1)

samp_accessions="${base_dir}/analysis/input_data/pfalciparum/pacb_ilmn/data_accessions.tsv"

if [[ ! -e $samp_accessions ]]; then echo "$samp_accessions" required but not found; exit 1; fi

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
    dl_command="enaDataGet -f fastq -d ${output_dir}/${samp_name} $ilmn_read_accession"
    ## run as single job 
    eval $dl_command

    ## run on LSF cluster
    #submission_command="bsub -R \"select[mem>1000] rusage[mem=1000]\" -M1000 -o ${log_dir}/${samp_name}_ilmn_reads.o -e ${log_dir}/${samp_name}_ilmn_reads.e -J ${samp_name}_ena_dl "
    #eval $submission_command $dl_command
done

