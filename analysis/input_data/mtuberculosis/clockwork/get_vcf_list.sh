[[ $# != 1 ]] && echo "usage: $0 num_vcfs" && exit 1
num_vcfs=$1
clockwork_vcf_list="clockwork_vcf_list.txt"

if [[ -e $clockwork_vcf_list ]]; then echo "$clockwork_vcf_list already exists, skipping its production"
else
	for base_dir in $(cut -f10 clockwork_mykrobe_samples/clockwork_vcfs_dirs.tsv)
	do
		find ${base_dir}/cortex -name "*FINAL*raw.vcf" >> $clockwork_vcf_list
	done
fi

indexed_vcf_dir="/hps/nobackup2/zi/bletcher/vcfs/tb/clockwork_vcfs"
mkdir -p $indexed_vcf_dir

indexed_vcf_list="indexed_vcf_list.txt"
> $indexed_vcf_list
i=0
for vcf in $(cat $clockwork_vcf_list); do
	i=$((i + 1))
	fname=${indexed_vcf_dir}/s${i}.vcf
	cp $vcf $fname
	echo ${fname}.gz >> $indexed_vcf_list
	bgzip $fname
	bcftools index ${fname}.gz
done


cat <(head -n $num_vcfs $indexed_vcf_list) vcfs_Comas_pacb_ilmn.txt > "vcf_list.txt"


