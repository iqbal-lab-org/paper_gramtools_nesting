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
    metric="Average Memory"
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
mean_reads=$(get_metric "all reads:" "tb_gramtools_genotype")

get_speed(){
    subdir1=$1
    subdir2=$2
    metric="CPU time"
    time1=$(get_metric "$metric" "$subdir1")
    time2=$(get_metric "$metric" "$subdir2")
    total_time=$(echo "print($mean_reads / ($time1 + $time2))" | python3)
    echo "$total_time"
}


### Genotyping ###

# Get mean CPU times for genotyping
mean_CPU_gramtools=$(get_metric "CPU time" "tb_gramtools_genotype")
mean_CPU_vg=$(get_speed tb_vg_map tb_vg_genotype)
mean_CPU_graphtyper=$(get_speed tb_map_reads_to_ref tb_graphtyper_genotype)


# Get mean RAM use
mean_RAM_gramtools=$(get_RAM "tb_gramtools_genotype")
mean_RAM_vg=$(get_RAM "tb_vg_map" "tb_vg_genotype")
mean_RAM_graphtyper=$(get_RAM "tb_map_reads_to_ref" "tb_graphtyper_genotype")

##### Build/indexing (vg and gramtools)
mean_RAM_build_gramtools=$(get_RAM "tb_gram_build")
mean_CPU_build_gramtools=$(get_metric "CPU time" "tb_gram_build")
mean_RAM_build_vg=$(get_RAM "tb_vg_build")
mean_CPU_build_vg=$(get_metric "CPU time" "tb_vg_build")

echo -e "Tool\tOperation\tAverage RAM\tSpeed(reads/CPU sec or CPU sec)" > $ofile
echo -e "vg\tindex\t$mean_RAM_build_vg\t$mean_CPU_build_vg" >> $ofile
echo -e "gramtools\tindex\t$mean_RAM_build_gramtools\t$mean_CPU_build_gramtools" >> $ofile
echo -e "vg\tgenotype\t$mean_RAM_vg\t$mean_CPU_vg" >> $ofile
echo -e "gramtools\tgenotype\t$mean_RAM_gramtools\t$mean_CPU_gramtools" >> $ofile
echo -e "graphtyper\tgenotype\t$mean_RAM_graphtyper\t$mean_CPU_graphtyper" >> $ofile
