import json

import click

from jvcf_processing import Region, get_region, evaluate_site, get_sites_in_region

result_fields = [
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
    "ambiguous",
    "is_parent_site",
]
result_template = {k: "NA" for k in result_fields}


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
    callback=get_region,
)
def main(genotyped_jvcf, truth_jvcf, region: Region, output_file):
    """
    :genotyped_jvcf: A single sample jvcf which either has same sites as :truth_jvcf: or has sites in :region: corresponding to :truth_jvcf:
    """
    with open(genotyped_jvcf) as fone, open(truth_jvcf) as ftwo:
        genotyped = json.load(fone)
        truth = json.load(ftwo)

    truth_sites = truth["Sites"]
    genotyped_sites = genotyped["Sites"]
    genotyped_sites = get_sites_in_region(genotyped_sites, region)
    if len(genotyped_sites) != len(truth_sites):
        raise ValueError(
            f"{len(genotyped_sites)} genotyped sites vs {len(truth_sites)} truth sites, should be same number. Use --region ?"
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

        if str(i) in truth["Child_Map"]:
            next_result["is_parent_site"] = 1
        else:
            next_result["is_parent_site"] = 0

        # Make sure no new keys will be introduced
        next_result.update(
            {key: val for key, val in eval_results.items() if key in next_result}
        )
        print("\t".join(map(str, next_result.values())) + "\n")
        # fout.write("\t".join(map(str, next_result.values())) + "\n")


if __name__ == "__main__":
    main()
