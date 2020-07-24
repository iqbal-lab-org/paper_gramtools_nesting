for vcf in $(find ./analysis/outputs/pacb_ilmn_validation/pf_4surfants/465102ee/gramtools_genotype -name "*.vcf.gz"); do
	nozip=${vcf%*.gz}
	bgzip -d "$vcf"
	sed -i '/ID=COV/s/Type=Integer/Type=Float/' "$nozip"
	bgzip "$nozip"
done
