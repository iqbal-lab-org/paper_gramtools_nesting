usage(){
	echo "usage: $0 input.tsv visualise_prg_script coverage_graph output_dir"
	exit 1
}

input_tsv=$1
vis_executable=$2
cov_graph=$3
output_dir=$4

if [[ -z $input_tsv || -z $vis_executable || -z $cov_graph || -z $output_dir ]]; then usage; fi

mkdir -p $output_dir

IFS=$'\n'
for line in $(cat $input_tsv | tail -n+2); do
	IFS=$'\t'; elems=($line)
	if [[ ${elems[4]} == 0 || ${elems[3]} == " " ]]; then continue; fi
	ofname="${output_dir}/dot_site_${elems[0]}"
	$vis_executable $cov_graph ${elems[0]} ${elems[0]} $ofname
	cat "${ofname}.gv" | dot -Tpng > "${ofname}.png"
done
	
