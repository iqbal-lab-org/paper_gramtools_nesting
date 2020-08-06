set -eu
if [[ $# -ne 3 ]]; then echo "usage: $0 vcf_fofn min_del_size ofile"; exit 1; fi
vcf_fofn=$1
files=$(cat $vcf_fofn | tr '\n' ' ')
min_del_size=$2
ofile=$3

awk '{if (length($4) > '$min_del_size' && length($5) == 1) {print $0}}' $files | grep '1/1' | grep 'PASS'| awk '{print $1"\t"$2"\t"$2 + length($4)}' | bedtools sort | bedtools merge | awk '{l++;print $0"\tbigdel_"l}' > $ofile
