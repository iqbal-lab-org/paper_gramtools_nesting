import sys
import json

from jvcf_processing import jVCF_to_VCF


def usage():
    print(f"Usage: {sys.argv[0]} input_jvcf output_vcf_fname")
    exit(0)


if len(sys.argv) != 3:
    usage()

with open(sys.argv[1]) as fin, open(sys.argv[2], "w") as fout:
    input_jvcf = json.load(fin)
    converter = jVCF_to_VCF()
    converter.convert(input_jvcf, fout)
