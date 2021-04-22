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


def is_symbolic(allele: str) -> bool:
    return allele.startswith("<DEL") or allele.startswith("<INS")


def symbolic_del_to_sequence(vcf_record, ref_sequences):
    assert len(vcf_record.alleles) == 2
    pos = vcf_record.pos - 1
    event_size = vcf_record.info["SVLEN"]
    # Extracts the first base before the deletion + the deleted sequence
    ref_seq = ref_sequences[vcf_record.chrom][pos - 1 : pos + event_size]
    alt_seq = ref_seq[0]
    vcf_record.alleles = [ref_seq, alt_seq]
    vcf_record.pos = pos
    return vcf_record


def symbolic_ins_to_sequence(vcf_record, ref_sequences):
    assert len(vcf_record.alleles) == 2
    pos = vcf_record.pos - 1
    event_size = vcf_record.info["SVLEN"]
    inserted_seq = vcf_record.info["SEQ"]
    ref_seq = ref_sequences[vcf_record.chrom][pos]
    alt_seq = ref_seq + inserted_seq
    vcf_record.alleles = [ref_seq, alt_seq]
    return vcf_record


def process_symbolic(record, ref_sequences):
    first_alt = record.alts[0]
    assert is_symbolic(first_alt)
    if first_alt.startswith("<DEL"):
        return symbolic_del_to_sequence(record, ref_sequences)
    else:
        return symbolic_ins_to_sequence(record, ref_sequences)


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
@click.option("--remove_ref_and_null", is_flag=True)
def main(fasta_ref, vcf_file, output_vcf, remove_ref_and_null):
    output_fname = Path(output_vcf)
    output_fname.parent.mkdir(exist_ok=True, parents=True)

    ref_seqs = {}
    for seq in pyfastaq.sequences.file_reader(fasta_ref):
        seq_id = seq.id.split()[0]
        ref_seqs[seq_id] = seq.seq

    with open(vcf_file) as fin, open(output_vcf, "w") as fout:
        varfile = VariantFile(fin)
        fout.write(str(varfile.header))
        found_symbolics = set()
        for record in varfile:
            if remove_ref_and_null:
                # Ignore ref calls; these force varifier to unnecessarily ignore overlaps
                first_gt_call = record.samples[0]["GT"][0]
                if first_gt_call == 0 or first_gt_call is None:
                    continue
            first_alt = record.alts[0]
            if is_symbolic(first_alt):
                if "AGGREGATED" in first_alt:
                    assert record.id not in found_symbolics
                    found_symbolics.add(record.id)
                    record = process_symbolic(record, ref_seqs)
                else:
                    assert record.id.rsplit(".", maxsplit=1)[0] in found_symbolics
                    continue
            fout.write(str(record))
        print(
            f"Processed {len(found_symbolics)} symbolically encoded structural variants"
        )


if __name__ == "__main__":
    main()
