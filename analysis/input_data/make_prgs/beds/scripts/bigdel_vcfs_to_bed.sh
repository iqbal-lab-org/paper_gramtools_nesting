set -eu
if [[ $# -ne 3 ]]; then echo "usage: $0 vcf_fofn min_del_size ofile"; exit 1; fi
vcf_fofn=$1
files=$(cat $vcf_fofn | tr '\n' ' ')
min_del_size=$2
ofile=$3

origin_ofile=${ofile}.origin
echo -e "contig\tstart\tstop\tfound_in_samples" > $origin_ofile


for fname in $files; do
	gunzip -c $fname | awk '{if (length($4) > '$min_del_size' && length($5) == 1) {print $0}}' | grep '1/1' | grep 'PASS'| awk '{print $1"\t"$2 - 1"\t"$2 + length($4) - 1"\t'$(basename ${fname%%.*})'"}'
	done | bedtools sort | tee -a $origin_ofile | bedtools merge -o distinct -c 4 | cut -f1-3 | awk '{l++;print $0"\tbigdel_"l}' > $ofile
