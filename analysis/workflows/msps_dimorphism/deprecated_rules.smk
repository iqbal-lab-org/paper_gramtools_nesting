rule tree_basal_split:
    input:
        f"{output_trees}/{{gene}}/RAxML_bestTree.{{gene}}",
    output:
        split=f"{output_tree_splits}/{{gene}}.basal_split",
        plot=f"{output_tree_splits}/{{gene}}.pdf",
    params:
        script=f'{config["scripts"]}/{WORKFLOW}/tree_partition.py',
        output_prefix=f"{output_tree_splits}/{{gene}}",
    shell:
        "python3 {params.script} {input} {params.output_prefix}"


rule build_graphs:
    input:
        res_json=f"{output_gtyping}/combined.json",
        genes_bed=config["genes_bed"],
        tree_split=f"{output_tree_splits}/{{gene}}.basal_split",
    output:
        all_samples=f"{output_graphs}/{{gene}}_all_samples.gml",
        partition1=f"{output_graphs}/{{gene}}_partition_1.gml",
        partition2=f"{output_graphs}/{{gene}}_partition_2.gml",
    params:
        output_prefix=f"{output_graphs}/{{gene}}",
        script=f'{config["scripts"]}/{WORKFLOW}/get_site_diversity_graphs.py',
    shell:
        """
        # Produce region
        match=$(grep -w {wildcards.gene} {input.genes_bed})
        if [[ -z ${{match}} ]]; then echo "ERROR: No match of gene name in bed file"; exit 1; fi
        IFS="\t"; elems=($match)    
        adj_start=$((${{elems[1]}} + 1))
        reg="${{elems[0]}}:${{adj_start}}-${{elems[2]}}"

        python3 {params.script} {input.res_json} {params.output_prefix} --region $reg # for all_samples
        python3 {params.script} {input.res_json} {params.output_prefix} --region $reg -p {input.tree_split} # for partitions
        """
