rule minimap2_call:
	input:
		assembly=get_assembly,
		ref=config["starting_prg"]["fasta_ref"]
	output:
		gzipped=f'{output_minimap2}/{{sample}}_vars.vcf.gz',
		indexed=f'{output_minimap2}/{{sample}}_vars.vcf.gz.csi',
	params:
		vcf=f'{output_minimap2}/{{sample}}_vars.vcf',
	shell:
		"""
		minimap2 -c --cs {input.ref} {input.assembly} | sort -k6,6n -k8,8n | paftools.js call -l50 -L50 -f {input.ref} -s {wildcards.sample} - > {params.vcf}
		bgzip {params.vcf} && bcftools index {output.gzipped}
		"""

rule validate_input_regions:
	input:
		vcfs=expand(f'{output_minimap2}/{{sample}}_vars.vcf.gz',sample=SAMPLES),
		input_regions=config["beds"]["with_origin"],
	output:
		vcf=f'{output_regions}/minimap2_vars.vcf.gz',
		tsv=f'{output_regions}/minimap2_validated.tsv',
	params:
		validation_script=f'{config["scripts"]}/{WORKFLOW}/find_input_dels.py'
	shell:
		"""
		bcftools merge {input.vcfs} -Oz -o {output.vcf}
		python3 {params.validation_script} {output.vcf} {input.input_regions} {output.tsv}
		"""


rule merge_and_assess_vcfs:
	input:
		vcfs=expand(f'{output_genotyped}/{{condition}}/{{sample}}.vcf.gz', sample = SAMPLES, allow_missing=True),
		input_regions=config["beds"]["with_origin"],
	output:
		vcf=f'{output_genotyped}/{{condition}}_merged.vcf.gz',
		tsv=f'{output_genotyped}/{{condition}}_validated.tsv',
	params:
		validation_script=f'{config["scripts"]}/{WORKFLOW}/find_input_dels.py'
	shell:
		"bcftools merge {input.vcfs} -i - -m all -Oz -o {output.vcf};"
		"python3 {params.validation_script} {output.vcf} {input.input_regions} {output.tsv}"
