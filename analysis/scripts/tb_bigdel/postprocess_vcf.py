"""
VCF postprocessing:
* Convert symbolic SVs (<DEL>, <INS>) into ref/alt sequence. These occur in graphtyper2
  output VCFs and cannot be analysed by varifier (nor by bcftools consensus for <INS>)
* Optionally remove REF and null calls: if left in these can cause varifier to unnecessarily 
  ignore some overlapping variant records
"""
import sys
from pathlib import Path

import click
import pyfastaq
from pysam import VariantFile

@click.command()
@click.argument(
    "fasta_ref",
    type=click.Path(exists=True),
)
@click.argument(
    "vcf_file",
    type=click.Path(exists=True),
)
@click.argument("output_vcf")
@click.option("--remove_ref_and_null",is_flag=True)
def main(fasta_ref, vcf_file, output_vcf, remove_ref_and_null):
    output_fname = Path(output_vcf)
    output_fname.parent.mkdir(exist_ok=True, parents=True)

    ref_seqs = {}
    for seq in pyfastaq.sequences.file_reader(fasta_ref):
        seq_id = seq.id.split()[0]
        ref_seqs[seq_id] = seq.seq

    with open(vcf_file) as fin, open(output_vcf,"w") as fout:
        varfile = VariantFile(fin)
        fout.write(str(varfile.header))
        for record in varfile:
            if remove_ref_and_null:
                # Ignore ref calls; these force varifier to unnecessarily ignore overlaps
                first_gt_call = record.samples[0]["GT"][0]
                if first_gt_call == 0 or first_gt_call is None:
                    continue
            #first_alt = record.alts[0]
            #if first_alt.startswith("<DEL"):
            #    del_size = record.info["SVLEN"]
            fout.write(str(record))
if __name__ == "__main__":
    main()
