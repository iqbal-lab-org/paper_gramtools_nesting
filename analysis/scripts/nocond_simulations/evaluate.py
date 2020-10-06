import json
import sys
from typing import NamedTuple, Union

import pandas as pd
import click

from jvcf_processing import (
    find_sample_index,
    AlleleCall,
    get_called_allele,
    evaluate_site,
    num_sites_under,
)

columns = [
    "prg",
    "simu_path",
    "err_rate",
    "fcov",
    "nesting",
    "res_has_call",
    "truth_has_call",
    "res_is_correct",
    "classif",
    "lvl_1",
    "GC",
    "GCP",
    "gt_allele",
    "truth_allele",
    "edit_dist",
    "cov_gt_allele",
    "cov_other_alleles",
    "site_num",
    "site_pos",
    "genotyped_ambiguous",
    "truth_ambiguous",
    "num_child_sites",
]
result_template = {k: "NA" for k in columns}


def print_cols(ctx, param, value):
    if not value:
        return
    print("\t".join(columns))
    ctx.exit()


@click.command()
@click.option(
    "-p",
    help="print tsv columns and exit",
    is_flag=True,
    callback=print_cols,
    expose_value=False,
    is_eager=True,
)
@click.option("-n", "--prg_name", required=True, type=str)
@click.option("--num", help="simu path number", required=True, type=int)
@click.option("-e", "--err_rate", required=True, type=int)
@click.option("-c", "--fcov", required=True, type=int)
@click.option("--nesting", default="", type=str)
@click.argument("truth_json", type=click.Path(exists=True))
@click.argument("res_json", type=click.Path(exists=True))
@click.argument("output_path")
def main(prg_name, num, err_rate, fcov, nesting, truth_json, res_json, output_path):
    result_template["prg"] = prg_name
    result_template["simu_path"] = str(num)
    result_template["err_rate"] = str(err_rate)
    result_template["fcov"] = str(fcov)
    result_template["nesting"] = nesting

    ## Load up truth json
    with open(truth_json) as fin:
        truth_json = json.load(fin)

    ## Load up result json
    with open(res_json) as fin:
        res_json = json.load(fin)
        lvl1_sites = set(res_json["Lvl1_Sites"])

    ## Evaluate calls
    fout = open(output_path, "w")
    # Below sample_id assumes gramtools simulate was called with that --sample_id
    sample_id = f"{prg_name}{num}"
    truth_sample_index = find_sample_index(truth_json, sample_id)
    assert sample_id == res_json["Samples"][0]["Name"]

    for i in range(len(truth_json["Sites"])):
        next_result = result_template.copy()

        called_site_json = res_json["Sites"][i]
        truth_site_json = truth_json["Sites"][i]

        eval_results = evaluate_site(
            called_site_json, 0, truth_site_json, truth_sample_index
        )
        # Make sure no new keys will be introduced
        next_result.update(
            {key: val for key, val in eval_results.items() if key in next_result}
        )

        next_result["num_child_sites"] = num_sites_under(res_json["Child_Map"], str(i))

        if lvl1_sites == {"all"} or i in lvl1_sites:
            next_result["lvl_1"] = "1"
        else:
            next_result["lvl_1"] = "0"

        next_result["site_num"] = i
        next_result["site_pos"] = called_site_json["POS"]

        fout.write("\t".join(map(str, next_result.values())) + "\n")

    fout.close()


if __name__ == "__main__":
    main()
