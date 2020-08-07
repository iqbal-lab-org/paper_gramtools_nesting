[[ $# != 1 ]] && echo "usage: $0 num_vcfs" && exit 1
num_vcfs=$1
vcf_list="clockwork_vcf_list.txt"
> $vcf_list
for base_dir in $(cut -f10 clockwork_mykrobe_samples/clockwork_vcfs_dirs.tsv)
do
	find ${base_dir}/cortex -name "*FINAL*raw.vcf" >> $vcf_list
done

cat <(head -n $num_vcfs $vcf_list) vcfs_Comas_pacb_ilmn.txt > "vcf_list.txt"

