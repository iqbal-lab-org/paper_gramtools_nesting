####
# Script to find computational performance of vg, graphtyper2 and gramtools on tuberculosis benchmark.
# To obtain disk use of 'index' operation, inspect files sizes in analysis/outputs/vg_builds and analysis/outputs/gram_builds
# Assumptions:
#   - log directories are named after rules
#   - patterns output by cluster scheduler in log files: eg 'CPU time' for CPU runtime
# These assumptions only hold if running on a LSF cluster system and with the snakemake lsf-profile.
 
set -e 

usage(){
    echo "usage: $0 log_dir ofile"
    exit 0
}

if [[ $# != 2 ]]; then usage; fi

log_dir="$1"
ofile="$2"

if [[ ! -e "$log_dir" ]]; then usage; fi

get_metric(){
    metric=$1
    log_sub_dir=$2
    found=$(find "${log_dir}/${log_sub_dir}" -name "*.out" | xargs grep "$1" | grep -Eo " [0-9]+" | awk '{tot+=$1}END{printf "%f", tot/NR}')
    echo $found
}

get_RAM(){
    metric="Max Memory"
    subdir1=$1
    RAM1=$(get_metric "$metric" "$subdir1")
    printf -v RAM1 %.0f "$RAM1"
    if [[ $# -eq 1 ]]; then
        echo "$RAM1"
    else
        subdir2=$2
        RAM2=$(get_metric "$metric" "$subdir2")
        printf -v RAM2 %.0f "$RAM2"
    echo "${RAM1}(${RAM2})"
    fi
}

# Get mean number of mapped reads
mean_reads=$(get_metric "all reads:" "msps_genotype")


### Genotyping ###

# Get mean CPU times for genotyping
mean_CPU_gramtools=$(get_metric "CPU time" "msps_genotype")
mean_CPU_gramtools=$(echo "print($mean_reads / $mean_CPU_gramtools)" | python3)


# Get mean RAM use
mean_RAM_gramtools=$(get_RAM "msps_genotype")

echo -e "Tool\tOperation\tMax RAM\tSpeed(reads/CPU sec or CPU sec)" > $ofile
echo -e "gramtools\tgenotype\t$mean_RAM_gramtools\t$mean_CPU_gramtools" >> $ofile
