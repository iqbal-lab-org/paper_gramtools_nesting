rule gram_build:
	input:
		prg=config["starting_prg"]["gram_prg"],
		ref=config["starting_prg"]["fasta_ref"]
		
	output:
		f'{output_gram_build}/cov_graph',

	params:
		k=config["starting_prg"]["gram_kmer_size"],
		gram_dir=output_gram_build,

	resources:
		mem_mb=20000

	shell:
		"""
		gramtools build --prg {input.prg} --ref {input.ref} --kmer_size {params.k} --gram_dir {params.gram_dir} --force
		"""

rule gramtools_genotype:
	input:
		gram_build_completed=rules.gram_build.output,
		gram_dir=f'{Path(rules.gram_build.output[0]).parent}',
		reads_files=get_reads

	output:
		sample_geno_dir=directory(f'{output_genotyped}/{conditions[0]}/{{sample}}'),
		gzipped=f'{output_genotyped}/{conditions[0]}/{{sample}}.vcf.gz',
		indexed=f'{output_genotyped}/{conditions[0]}/{{sample}}.vcf.gz.csi',
		jvcf=f'{output_genotyped}/{conditions[0]}/{{sample}}/genotype/genotyped.json',

	threads: 10
	resources:
		mem_mb=10000

	shell:
		"""
		gramtools genotype -i {input.gram_dir} -o {output.sample_geno_dir} --reads {input.reads_files} --sample_id {wildcards.sample} --max_threads {threads} --force
		cp {output.sample_geno_dir}/genotype/genotyped.vcf.gz {output.gzipped}
		bcftools index {output.gzipped}
		bcftools filter -i 'FT!="AMBIG"' {output.gzipped} -Oz -o tmp.vcf.gz
		mv tmp.vcf.gz {output.gzipped} && bcftools index -f {output.gzipped}
		"""


rule vg_build:
	input:
		vg_graph=config["starting_prg"]["vg_prg"]
	output:
		xg=f'{output_vg_build / "prg.xg"}',
		gcsa=f'{output_vg_build / "prg.gcsa"}',
		snarls=f'{output_vg_build / "prg.snarls"}',
	params:
		k=config["starting_prg"]["vg_kmer_size"],
		X=config["starting_prg"]["vg_doubling_steps"],
	shadow:
		"shallow"
	resources:
		mem_mb=25000
	threads:
		10
	shell:
		"""
		vg index -k {params.k} -x {output.xg} -L {input.vg_graph} 
		vg prune -r -t {threads} {input.vg_graph} > pruned.vg
		vg index -L -X {params.X} -k {params.k} -g {output.gcsa} -p -t {threads} pruned.vg
		vg snarls {output.xg} > {output.snarls}
		#vg deconstruct -p ref -e {input.vg_graph} {input.vg_graph} > {{output.vcf}}
		"""

rule vg_map:
	input:
		xg=rules.vg_build.output.xg,
		gcsa=rules.vg_build.output.gcsa,
		reads_files=get_reads,
	output:
		mapped_packed=f'{output_vg_mapped}/mapped_{{sample}}.pack',
	params:
		mapped_gam=f'{output_vg_mapped}/mapped_{{sample}}.gam',
	threads: 10
	resources:
		mem_mb=10000
	shadow:
		"shallow"
	shell:
		"""
		rfilecmd=""
		for rfile in {input.reads_files}; do rfilecmd="$rfilecmd -f $rfile"; done
		vg map -x {input.xg} -g {input.gcsa} $rfilecmd -t {threads} > {params.mapped_gam}
		vg pack -x {input.xg} -g {params.mapped_gam} -Q 5 -o {output.mapped_packed}
		"""

rule vg_genotype:
	input:
		xg=rules.vg_build.output.xg,
		mapped=rules.vg_map.output.mapped_packed,
		snarls=rules.vg_build.output.snarls,
		vcf_to_genotype=config["vcf_to_genotype"]
	output:
		gzipped=f'{output_genotyped}/{conditions[1]}/{{sample}}.vcf.gz',
		indexed=f'{output_genotyped}/{conditions[1]}/{{sample}}.vcf.gz.csi',
	params:
		vcf=f'{output_genotyped}/{conditions[1]}/{{sample}}.vcf',
	shell:
		"""
		vg call {input.xg} -k {input.mapped} -v {input.vcf_to_genotype} -r {input.snarls} -s {wildcards.sample} --ploidy 1 > {params.vcf}
		bgzip {params.vcf} && bcftools index {output.gzipped}
		"""

rule baseline_ref_genotype:
	"""Dummy rule producing empty vcfs used for baseline gene portion-assembly alignments"""
	input:
		vcf=config["vcf_template"]
	output:
		gzipped=temp(f'{output_genotyped}/{conditions[2]}/{{sample}}.vcf.gz'),
		indexed=f'{output_genotyped}/{conditions[2]}/{{sample}}.vcf.gz.csi',
	params:
		vcf=f'{output_genotyped}/{conditions[2]}/{{sample}}.vcf',
	shell:
		"""
		cp {input.vcf} {params.vcf}
		sed -i 's/sample/{wildcards.sample}/' {params.vcf}
		bgzip {params.vcf} && bcftools index {output.gzipped}
		"""


rule graphtyper_genotype:
	input:
		bam=f'{output_ref_alignments}/{{sample}}.bam',
		var_regions=config["beds"]["with_flank"],
		vcf=config["vcf_to_genotype"],
		ref=config["starting_prg"]["fasta_ref"],
	output:
		gzipped=f'{output_genotyped}/{conditions[3]}/{{sample}}.vcf.gz',
		indexed=f'{output_genotyped}/{conditions[3]}/{{sample}}.vcf.gz.csi',
	shadow:
		"shallow"
	threads: 10
	resources:
		mem_mb=10000
	shell:
		"""
		output_dir=tmp_{wildcards.sample}
		mkdir -p $output_dir
		region_file=${{output_dir}}/regions.txt
		> $region_file
		IFS="\n"; for gene_line in $(cat {input.var_regions})
		do
			ref_name=${{elems[0]}}
			IFS="\t"; elems=($gene_line)	
			adjusted_start=$((${{elems[1]}} + 1))
			reg="${{elems[0]}}:${{adjusted_start}}-${{elems[2]}}"
			echo $reg >> $region_file
		done
		graphtyper genotype_sv --vverbose --region_file $region_file --output $output_dir --sam {input.bam} --threads {threads} {input.ref} {input.vcf}

		# Concat in positionally sorted order, and make DEL and INS symbolic alleles applyable by bcftools consensus downstream
		IFS=$'\n'
		vcfs_made=($(find ${output_dir}/${ref_name} -name "*.vcf.gz" | sort -n))
		bcftools concat "${vcfs_made[@]}" | sed 's/<DEL.*>/<DEL>/' | sed 's/<INS.*>/<INS>/' > tmp.vcf

		# Change sample name
		echo "$(bcftools query -l tmp.vcf) {wildcards.sample}" > rename.txt
		bcftools reheader -s rename.txt tmp.vcf > renamed.vcf
		bgzip -c renamed.vcf > {output.gzipped}
		bcftools index {output.gzipped}
		"""