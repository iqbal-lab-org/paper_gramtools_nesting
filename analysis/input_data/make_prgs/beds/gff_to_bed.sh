set -eu
genes_file=$1
gff_file=$2
output_bed=$3

[[ $# != 3 ]] && echo "usage: $0 genes_file gff_file output_bed" \
	&& echo -e "genes_file: one gene name per line" \
	&& echo -e "gff_file: genome annotation file" \
	&& exit 0

genes=($(cat $genes_file))
grep_cmd="grep $gff_file"
> $output_bed
for gene in ${genes[@]}
	do 
		found_annot=$($grep_cmd -e "$gene;")
		if [[ $(echo $found_annot | wc -l) -gt 1 ]]; then
			echo "Error: found more than one annotation line for gene $gene"
			exit 1
		fi
		echo $found_annot | awk '{$4-=1;print $1"\t"$4"\t"$5"\t""'$gene'"}' >> $output_bed
done
