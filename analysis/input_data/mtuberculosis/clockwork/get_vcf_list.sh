### From a file describing directories containing
### vcf files made by cortex(and clockwork) on tuberculosis,
### copies them locally, indexes them and lists them

set -eu

[[ $# != 2 ]] && echo "usage: $0 num_vcfs base_directory (the root of this repository)" && exit 1
num_vcfs=$1
base_dir=$(realpath $2)
base_dir="${base_dir}/analysis/input_data/mtuberculosis/clockwork"

clockwork_vcf_list="${base_dir}/intermediate_file_lists/clockwork_vcf_list.txt"

if [[ -e $clockwork_vcf_list ]]; then echo "$clockwork_vcf_list already exists, skipping its production"
else
	for vcf_base_dir in $(cut -f10 ${base_dir}/clockwork_mykrobe_samples/clockwork_vcfs_dirs.tsv)
	do
		find ${vcf_base_dir}/cortex -name "*FINAL*raw.vcf" >> $clockwork_vcf_list
	done
fi

indexed_vcf_dir="${base_dir}/vcfs"
mkdir -p $indexed_vcf_dir

indexed_vcf_list="${base_dir}/intermediate_file_lists/indexed_vcf_list.txt"
> $indexed_vcf_list
i=0
for vcf in $(cat $clockwork_vcf_list); do
    i=$((i + 1))
    if [[ $i -gt $num_vcfs ]];then break; fi
    fname=${indexed_vcf_dir}/s${i}.vcf
    gzipped_fname="${fname}.gz"
    if [[ ! -e "$gzipped_fname" ]]; then
        cp $vcf $fname
        bgzip $fname
        bcftools index $gzipped_fname
    fi
    recorded_fname=$(basename $gzipped_fname)
    echo $recorded_fname >> $indexed_vcf_list
done


cat <(head -n $num_vcfs $indexed_vcf_list) ${base_dir}/intermediate_file_lists/vcfs_Comas_pacb_ilmn.txt > "${base_dir}/vcf_list.txt"


