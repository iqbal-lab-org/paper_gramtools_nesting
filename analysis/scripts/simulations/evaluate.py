import json
import sys
from typing import Tuple

import pandas as pd
import click
import edlib

columns = [
    "prg",
    "simu_path",
    "err_rate",
    "fcov",
    "is_null_gt",
    "correct",
    "lvl_1",
    "GC",
    "GCP",
    "edit_dist",
]
result_template = {k: "None" for k in columns}

def find_sample_index(truth_json, sample_id: str) -> int:
    sample_index = -1
    for i in range(len(truth_json["Samples"])):
        if truth_json["Samples"][i]["Name"] == sample_id:
            sample_index = i
            break
    assert sample_index > -1
    return sample_index


def get_called_allele(sample_index: int, site_json) -> Tuple[str,str]:
    allele_index = site_json["GT"][sample_index]
    assert len(allele_index) == 1
    allele_index = allele_index[0]

    if allele_index == None:
        return "", "1"

    return site_json["ALS"][allele_index], "0"


def print_cols(ctx, param, value):
    if not value:
        return
    print("\t".join(columns))
    ctx.exit()

@click.command()
@click.option('-p',help='print tsv columns and exit',is_flag=True, callback=print_cols, expose_value=False, is_eager=True)
@click.option('-n','--prg_name',required=True, type=str)
@click.option('--num',help='simu path number',required=True, type=int)
@click.option('-e','--err_rate',required=True, type=int)
@click.option('-c','--fcov',required=True, type=int)
@click.argument('truth_json', type=click.Path(exists=True))
@click.argument('res_json', type=click.Path(exists=True))
@click.argument('output_path')
def main(prg_name, num, err_rate, fcov, truth_json, res_json, output_path):


    result_template["prg"] = prg_name
    result_template["simu_path"] = str(num)
    result_template["err_rate"] = str(err_rate)
    result_template["fcov"] = str(fcov)


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
    truth_sample_index = find_sample_index(truth_json, sample_id = f'{prg_name}{num}')

    for i in range(len(truth_json["Sites"])):
        next_result = result_template.copy()
        true_allele, true_null_call = get_called_allele(truth_sample_index, truth_json["Sites"][i])

        called_site_json = res_json["Sites"][i]
        called_allele, is_null_call = get_called_allele(0, called_site_json)

        next_result["is_null_gt"] = is_null_call

        if called_allele == true_allele and is_null_call == true_null_call:
            next_result["correct"] = "1"
        else:
            next_result["correct"] = "0"

        if lvl1_sites == {"all"} or i in lvl1_sites:
            next_result["lvl_1"] = "1"
        else:
            next_result["lvl_1"] = "0"

        next_result["GC"] = str(called_site_json["GT_CONF"][0])
        GCP_entry = called_site_json.get("GCP")
        if GCP_entry is not None:
            next_result["GCP"] = str(GCP_entry[0])

        next_result["edit_dist"] = str(
            edlib.align(called_allele, true_allele)["editDistance"]
        )

        fout.write("\t".join(next_result.values()) + "\n")

    fout.close()

if __name__ == "__main__":
    main()
