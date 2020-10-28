set -eu

[[ $# != 1 ]] && echo "usage: $0 base_directory (root path of the repository)" && exit 1

base_dir=$(realpath $1)

samp_accessions="${base_dir}/analysis/input_data/pfalciparum/pf3k/pf3k_release_5.tsv"
reads_dir="${base_dir}/analysis/input_data/pfalciparum/pf3k/fastqs"

#pf3k_ftp="ftp://ngs.sanger.ac.uk/production/pf3k/release_5/BAM/"
#pf3k_dl_script="${base_dir}/analysis/input_data/download_data/pf3k_dl_bams/pf3k_dl_one.sh"
pf3k_dl_script="${base_dir}/analysis/input_data/download_data/pf3k_dl_ilmn_one.sh"

log_dir="${base_dir}/analysis/logs/downloads"

mkdir -p "$log_dir"
mkdir -p "$reads_dir"

IFS=$'\n'
for line in $(cat ${samp_accessions})
do
	IFS=$'\t'; elems=($line)
    dl_this_sample="${elems[2]}"
    if [[ $dl_this_sample != '1' ]]; then
        continue
    else
        samp_name="${elems[0]}"
        ilmn_samp_accession="${elems[1]}"
        ## Command to download fastqs from ENA
        mkdir -p "${reads_dir}/${samp_name}"
        dl_command="bash $pf3k_dl_script $ilmn_samp_accession ${reads_dir}/${samp_name}"
        ## Command to download BAM from Sanger ftp
        #dl_command="bash $pf3k_dl_script ${pf3k_ftp}/${accession} $reads_dir"

        ## Run sequentially
        eval"$dl_command"

        ## Run in parallel, submitted to an LSF cluster
        #submission_command="bsub -R \"select[mem>1000] rusage[mem=1000]\" -M1000 -o ${log_dir}/${samp_name}_ilmn_reads.o -e ${log_dir}/${samp_name}_ilmn_reads.e -J ${samp_name}_ena_dl "
        #eval "$submission_command" "$dl_command"
    fi
done

