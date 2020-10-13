set -eu

[[ $# != 1 ]] && echo "usage: $0 base_directory (root path of the repository)" && exit 1
base_dir=$(realpath $1)
samp_accessions="${base_dir}/analysis/input_data/pfalciparum/pf3k/bam_list.txt"
bam_dir="${base_dir}/analysis/input_data/pfalciparum/pf3k/bams"

pf3k_ftp="ftp://ngs.sanger.ac.uk/production/pf3k/release_5/BAM/"

log_dir="${base_dir}/analysis/logs/downloads"

mkdir -p "$log_dir"
mkdir -p "$bam_dir"

for accession in $(cat ${samp_accessions})
do
    dl_command="wget -N ${pf3k_ftp}/${accession} -O ${bam_dir}/${accession}"

    ## run as single job 
    eval $dl_command

    ## submit to LSF cluster
    #submission_command="bsub -R \"select[mem>1000] rusage[mem=1000]\" -M1000 -o ${log_dir}/${accession}_ilmn_reads.o -e ${log_dir}/${accession}_ilmn_reads.e -J ${accession}_ena_dl "
    #eval $submission_command $dl_command
done

