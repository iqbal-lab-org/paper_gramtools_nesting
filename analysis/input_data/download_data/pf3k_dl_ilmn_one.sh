# Downloads fastqs from ena and subsamples reads files using seqtk if too large
set -eu

[[ $# != 2 ]] && echo "usage: $0 ena_sample_accession reads_dir" && exit 1

samp_accession=$1
reads_dir=$2

reads_files=($(find $reads_dir -name "*fastq.gz" -o -name "*fq.gz"))

if [[ "${#reads_files[@]}" -eq 0 ]]; then
    enaGroupGet -f fastq -d "$reads_dir" "$samp_accession"
    reads_files=($(find $reads_dir -name "*fastq.gz" -o -name "*fq.gz"))
fi

max_size=4500000000
max_size_cushioned=4000000000
for read_file in "${reads_files[@]}"
do
    fsize=$(gzip -dc "$read_file" | wc -c)
    if [[ $fsize -ge $max_size ]]; then
        sampling_ratio=$(echo "$max_size_cushioned / $fsize" | bc -l)
        tmp_name=${read_file%.gz}
        seqtk sample "$read_file" "$sampling_ratio" | gzip > "$tmp_name"
        mv "$tmp_name" "$read_file"
    fi
done
