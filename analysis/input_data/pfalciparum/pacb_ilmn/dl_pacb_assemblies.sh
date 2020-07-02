set -eu
base_dir="/nfs/leia/research/iqbal/bletcher/analyses/nesting_paper"

output_dir="${base_dir}/analysis/input_data/pfalciparum/pacb_ilmn/pacb_assemblies"
mkdir -p "$output_dir"

log_dir="${base_dir}/analysis/logs/downloads"
mkdir -p "$log_dir"

bsub.py 1 "dl_pacb_assemblies" -o "${log_dir}/dl_pacb_assemblies.o" -e "${log_dir}/dl_pacb_assemblies.e" wget "ftp://ftp.sanger.ac.uk/pub/project/pathogens/Plasmodium/falciparum/PF3K/ReferenceGenomes_Version1/GENOMES/*.fasta.gz" -nc -P "${output_dir}"
