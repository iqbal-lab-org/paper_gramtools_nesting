import json
import sys
from typing import NamedTuple, Union

import pandas as pd
import click
import edlib

columns = [
    "prg",
    "simu_path",
    "err_rate",
    "fcov",
    "nesting",
    "res_has_call",
    "truth_has_call",
    "res_is_correct",
    "lvl_1",
    "GC",
    "GCP",
    "edit_dist",
    "truth_allele",
    "gt_allele",
    "cov_gt_allele",
    "cov_other_alleles",
    "site_num",
    "site_pos"
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

class AlleleCall(NamedTuple):
    gt: Union[int, None]
    allele: str

    def has_call(self):
        return self.gt is not None


def get_called_allele(sample_index: int, site_json) -> AlleleCall: 
    allele_index = site_json["GT"][sample_index]
    assert len(allele_index) == 1
    allele_index = allele_index[0]

    allele = ""
    if allele_index is not None:
        allele = site_json["ALS"][allele_index]

    return AlleleCall(allele_index, allele)


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
@click.option('--nesting',default="", type=str)
@click.argument('truth_json', type=click.Path(exists=True))
@click.argument('res_json', type=click.Path(exists=True))
@click.argument('output_path')
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
    sample_id = f'{prg_name}{num}'
    truth_sample_index = find_sample_index(truth_json, sample_id)
    assert(sample_id == res_json["Samples"][0]["Name"])

    for i in range(len(truth_json["Sites"])):
        next_result = result_template.copy()
        true_call: AlleleCall = get_called_allele(truth_sample_index, truth_json["Sites"][i])

        called_site_json = res_json["Sites"][i]
        res_call = get_called_allele(0, called_site_json)

        next_result["res_has_call"] = res_call.has_call()
        next_result["truth_has_call"] = true_call.has_call()

        next_result["res_is_correct"] = res_call.allele == true_call.allele

        if lvl1_sites == {"all"} or i in lvl1_sites:
            next_result["lvl_1"] = "1"
        else:
            next_result["lvl_1"] = "0"

        next_result["GC"] = called_site_json["GT_CONF"][0]
        GCP_entry = called_site_json.get("GT_CONF_PERCENTILE")
        if GCP_entry is not None:
            next_result["GCP"] = GCP_entry[0]

        next_result["edit_dist"] = edlib.align(res_call.allele, true_call.allele)["editDistance"]
        next_result["truth_allele"] = true_call.allele
        next_result["gt_allele"] = res_call.allele
        if res_call.has_call():
            next_result["cov_gt_allele"] = called_site_json["COV"][0][res_call.gt]
            next_result["cov_other_alleles"] = called_site_json["DP"][0] - next_result["cov_gt_allele"]
        else:
            next_result["cov_other_alleles"] = called_site_json["DP"][0]
        next_result["site_num"] = i
        next_result["site_pos"] = called_site_json["POS"]

        fout.write("\t".join(map(str,next_result.values())) + "\n")

    fout.close()

if __name__ == "__main__":
    main()
