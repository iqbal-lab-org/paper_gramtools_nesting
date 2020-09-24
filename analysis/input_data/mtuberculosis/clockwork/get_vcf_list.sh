### From a file describing directories containing
### vcf files made by cortex(and clockwork) on tuberculosis,
### copies them locally, indexes them and lists them
###
### This is EBI cluster-specific

set -e 

[[ $# != 2 ]] && echo "usage: $0 num_vcfs cluster" && exit 1
num_vcfs=$1
cluster=$2

if [[ "$cluster" != "yoda" && "$cluster" != "noah" ]]; then echo "cluster should be yoda or noah"; exit 1; fi

clockwork_vcf_list="intermediate_file_lists/clockwork_vcf_list.txt"

if [[ -e $clockwork_vcf_list ]]; then echo "$clockwork_vcf_list already exists, skipping its production"
else
	for base_dir in $(cut -f10 clockwork_mykrobe_samples/clockwork_vcfs_dirs.tsv)
	do
		find ${base_dir}/cortex -name "*FINAL*raw.vcf" >> $clockwork_vcf_list
	done
fi

if [[ "$cluster" == "yoda" ]]; then
    indexed_vcf_dir="/hps/nobackup2/zi/bletcher/vcfs/tb/clockwork_vcfs"
else
    indexed_vcf_dir="/hps/nobackup/research/iqbal/bletcher/vcfs/tb/clockwork_vcfs"
fi
mkdir -p $indexed_vcf_dir

indexed_vcf_list="intermediate_file_lists/indexed_vcf_list.txt"
> $indexed_vcf_list
i=0
for vcf in $(cat $clockwork_vcf_list); do
    i=$((i + 1))
    fname=${indexed_vcf_dir}/s${i}.vcf
    gzipped_fname="${fname}.gz"
    if [[ ! -e "$gzipped_fname" ]]; then
        cp $vcf $fname
        bgzip $fname
        bcftools index $gzipped_fname
    fi
    echo $gzipped_fname >> $indexed_vcf_list
done


cat <(head -n $num_vcfs $indexed_vcf_list) intermediate_file_lists/vcfs_Comas_pacb_ilmn_"$cluster".txt > "vcf_list.txt"


