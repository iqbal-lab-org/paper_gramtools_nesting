set -eu

[[ $# != 1 ]] && echo "usage: $0 base_directory (root path of the repository)" && exit 1
base_dir=$(realpath $1)

output_dir="${base_dir}/analysis/input_data/pfalciparum/pacb_ilmn/pacb_assemblies"
mkdir -p "$output_dir"

log_dir="${base_dir}/analysis/logs/downloads"
mkdir -p "$log_dir"

dl_command="wget \"ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/*.fasta.gz\" -nc -P $output_dir"

## run as simple job
eval $dl_command

## run on LSF cluster
#submission_command="bsub -R \"select[mem>1000] rusage[mem=1000]\" -M1000 -o ${log_dir}/dl_pacb_assemblies.o -e ${log_dir}/dl_pacb_assemblies.e -J dl_pacb_assemblies "
#eval $submission_command $dl_command
