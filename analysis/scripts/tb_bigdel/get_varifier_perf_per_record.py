import click
from pathlib import Path
from operator import itemgetter
from collections import defaultdict

import edlib

from varifier.vcf_stats import format_dict_to_edit_dist_scores, _frs_from_vcf_record
from cluster_vcf_records import vcf_file_read


## Below function copied from varifier codebase
def per_record_stats_from_vcf_file(infile):
    """Gathers stats for each record in a VCF file.
    Returns a list of dictionaries of stats. One dict per VCF line.
    List is sorted by ref seq name (CHROM), then position (POS)"""
    stats = []
    wanted_keys = [
        "DP",
        "DPF",
        "FRS",
        "GT_CONF",
        "GT_CONF_PERCENTILE",
        "VFR_IN_MASK",
        "VFR_ED_RA",
        "VFR_ED_TR",
        "VFR_ED_TA",
        "VFR_ALLELE_LEN",
        "VFR_ALLELE_MATCH_COUNT",
        "VFR_ALLELE_MATCH_FRAC",
        "VFR_RESULT",
    ]
    key_types = {
        "DP": int,
        "DPF": float,
        "GT_CONF": float,
        "GT_CONF_PERCENTILE": float,
        "FRS": float,
        "VFR_IN_MASK": int,
        "VFR_ED_RA": int,
        "VFR_ED_TR": int,
        "VFR_ED_TA": int,
        "VFR_ALLELE_MATCH_FRAC": float,
        "VFR_ALLELE_LEN": int,
        "VFR_ALLELE_MATCH_COUNT": int,
    }
    header_lines, vcf_records = vcf_file_read.vcf_file_to_list(infile)
    for record in vcf_records:
        record_stats = {x: record.FORMAT.get(x, "NA") for x in wanted_keys}
        record_stats["FRS"] = _frs_from_vcf_record(record)
        record_stats["CHROM"] = record.CHROM
        record_stats["POS"] = record.POS + 1
        record_stats["ALS"] = [record.REF] + record.ALT
        record_stats["GT"] = record.FORMAT.get("GT", [None])
        for key, key_type in key_types.items():
            try:
                record_stats[key] = key_type(record_stats[key])
            except:
                pass

        stats.append(record_stats)

    stats.sort(key=itemgetter("CHROM", "POS"))
    return stats


def load_regions(bed_fname):
    result = defaultdict(list)
    if bed_fname is None:
        return result
    with open(bed_fname) as fin:
        for line in fin:
            entries = line.split("\t")
            result[entries[0]].append([int(entries[1]) + 1, int(entries[2])])
    return result


def is_in_regions(record_stats, regions):
    if len(regions) == 0:  # No filtering by regions
        return True
    chrom_regions = regions.get(record_stats["CHROM"])
    if chrom_regions is None:
        return False
    pos = int(record_stats["POS"])
    for region in chrom_regions:
        if region[0] <= pos <= region[1]:
            return True
    return False


def get_variant_type_and_size(record_stats):
    ref_allele = record_stats["ALS"][0]
    gtype_call = int(record_stats["GT"].split("/")[0])  # Assumes haploid
    called_allele = record_stats["ALS"][gtype_call]

    variant_size = edlib.align(ref_allele, called_allele, task="distance")[
        "editDistance"
    ]
    len_ref_allele, len_called_allele = len(ref_allele), len(called_allele)
    if len(ref_allele) == len(called_allele):
        variant_type = "SNP" if variant_size == 1 else "MNP"
    elif len(ref_allele) > len(called_allele):
        variant_type = "DEL"
    else:
        variant_type = "INS"

    return variant_type, variant_size


headers = [
    "POS",
    "CHROM",
    "Var_type",
    "Event_size",
    "Sample",
    "Tool",
    "Metric",
    "Classif",
    "Eddist_perf_num",
    "Eddist_perf_denum",
]


def print_headers(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    print("\t".join(headers))
    ctx.exit()


@click.command()
@click.argument(
    "input_dir",
    type=click.Path(exists=True),
)
@click.argument("output_tsv")
@click.option(
    "--header_only",
    is_flag=True,
    is_eager=True,
    expose_value=False,
    callback=print_headers,
    help="Produce the header for output_tsv and exit",
)
@click.option("--sample_name", type=str, required=True)
@click.option("--tool_name", type=str, required=True)
@click.option("--region_file", type=click.Path(exists=True), default=None)
def main(input_dir, output_tsv, sample_name, tool_name, region_file):

    precision_vcf = Path(input_dir) / "precision.vcf"
    recall_vcf = Path(input_dir) / "recall" / "recall.vcf"
    for fpath in [precision_vcf, recall_vcf]:
        assert fpath.exists()

    fout_path = Path(output_tsv)
    fout_path.parent.mkdir(exist_ok=True, parents=True)

    regions = load_regions(region_file)

    with fout_path.open("w") as fout:
        for vcf_fname, metric in zip(
            [precision_vcf, recall_vcf], ["precision", "recall"]
        ):
            all_stats = per_record_stats_from_vcf_file(str(vcf_fname))
            num_not_in_region = 0
            for record_stats in all_stats:
                if not is_in_regions(record_stats, regions):
                    num_not_in_region += 1
                    continue
                var_type, event_size = get_variant_type_and_size(record_stats)
                if event_size == 0:  # Can occur, eg AMBIG call
                    continue
                ed_num, ed_denum = format_dict_to_edit_dist_scores(record_stats)
                if ed_num is not None:
                    fout.write(
                        f'{record_stats["POS"]}\t{record_stats["CHROM"]}\t'
                        f"{var_type}\t{event_size}\t{sample_name}\t{tool_name}\t"
                        f'{metric}\t{record_stats["VFR_RESULT"]}\t'
                        f"{ed_num}\t{ed_denum}\n"
                    )
            print(f"{metric}:Skipped {num_not_in_region} variants not in {region_file}")


if __name__ == "__main__":
    main()
