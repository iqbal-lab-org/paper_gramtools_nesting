import json
import re
from pathlib import Path

import click

from jvcf_processing import (
    Region,
    click_get_region,
    evaluate_site,
    first_idx_in_region,
    get_n_sites_starting_from_region,
    num_sites_under,
)

result_fields = [
    "gene",
    "sample",
    "res_has_call",
    "truth_has_call",
    "res_is_correct",
    "classif",
    "GC",
    "GCP",
    "gt_allele",
    "truth_allele",
    "edit_dist",
    "cov_gt_allele",
    "cov_other_alleles",
    "site_num",
    "POS",
    "ambiguous",
    "num_child_sites",
]


def print_cols(ctx, param, value):
    if not value:
        return
    print("\t".join(result_fields))
    ctx.exit()


@click.command()
@click.option(
    "-p",
    help="print output tsv columns and exit",
    is_flag=True,
    callback=print_cols,
    expose_value=False,
    is_eager=True,
)
@click.argument("genotyped_jvcf", type=click.Path(exists=True))
@click.argument("truth_jvcf", type=click.Path(exists=True))
@click.argument("output_file", type=click.Path())
@click.option(
    "--region",
    "-r",
    help="In the form 'SEG:start-end', as in samtools/bcftools",
    default=None,
    callback=click_get_region,
)
def main(genotyped_jvcf, truth_jvcf, region: Region, output_file):
    """
    :genotyped_jvcf: A single sample jvcf which either has same sites as :truth_jvcf: or has sites in :region: corresponding to :truth_jvcf:
    """
    fout = open(output_file, "w")
    with open(genotyped_jvcf) as fone, open(truth_jvcf) as ftwo:
        genotyped = json.load(fone)
        truth = json.load(ftwo)

    truth_sites = truth["Sites"]
    genotyped_sites = genotyped["Sites"]
    site_num = first_idx_in_region(genotyped_sites, region)
    genotyped_sites = genotyped_sites[site_num : site_num + len(truth_sites)]
    if len(genotyped_sites) != len(truth_sites):
        raise ValueError(
            f"{len(genotyped_sites)} genotyped sites vs {len(truth_sites)} truth sites, should be same number. Use --region ?"
        )

    result_template = {k: "NA" for k in result_fields}
    fname_matcher = re.match("([^_]+)_([^_]+_[^_]+)_(.*).json", Path(truth_jvcf).name)
    result_template.update(
        {"sample": fname_matcher.groups()[0], "gene": fname_matcher.groups()[1]}
    )
    for i in range(len(truth_sites)):
        next_result = result_template.copy()

        called_site_json = genotyped_sites[i]
        truth_site_json = truth_sites[i]

        eval_results = evaluate_site(called_site_json, 0, truth_site_json, 0)
        # Make sure no new keys will be introduced
        next_result.update(
            {key: val for key, val in eval_results.items() if key in next_result}
        )

        next_result["POS"] = called_site_json["POS"]

        next_result["num_child_sites"] = num_sites_under(truth["Child_Map"], str(i))

        # Make sure no new keys will be introduced
        next_result.update(
            {key: val for key, val in eval_results.items() if key in next_result}
        )

        next_result["site_num"] = site_num
        site_num += 1
        # print("\t".join(map(str, next_result.values())) + "\n")
        fout.write("\t".join(map(str, next_result.values())) + "\n")
    fout.close()


if __name__ == "__main__":
    main()
