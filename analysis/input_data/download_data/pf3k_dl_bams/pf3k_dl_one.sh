# Downloads bam from ftp and subsamples it using samtools if too large
set -eu

[[ $# != 2 ]] && echo "usage: $0 ftp_location output_directory" && exit 1

ftp_location=$1
output_dir=$2
ofile=$(basename "$ftp_location")
ofile="$output_dir/$ofile"

#wget -N -c -P "${output_dir}" "${ftp_location}"
# Use this instead if the ftp server is not responsive enough and you trust you've already downloaded full file
if [[ ! -e "$ofile" ]]; then
    wget -N -c -P "${output_dir}" "${ftp_location}"
fi

fsize=$(wc -c <"$ofile")
max_size=4000000000
max_size_cushioned=3000000000
if [[ $fsize -ge $max_size ]]; then
    sampling_ratio=$(echo "$max_size_cushioned / $fsize" | bc -l | tr -d '.')
    samtools view -h -O BAM -o "${ofile}.tmp" -s 0."$sampling_ratio" "$ofile"
    mv "${ofile}.tmp" "$ofile"
fi
