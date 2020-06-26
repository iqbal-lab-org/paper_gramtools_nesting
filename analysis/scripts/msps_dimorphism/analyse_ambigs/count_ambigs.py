import sys
import json
import re


def usage():
    print(f"usage: {sys.argv[0]} jvcf.json output_fname [region]")
    exit(1)


if not 3 <= len(sys.argv) <= 4:
    usage()

jvcf_in = sys.argv[1]
with open(jvcf_in) as fin:
    jvcf = json.load(fin)

lvl1_sites = set(jvcf["Lvl1_Sites"])

output_fname = sys.argv[2]

region = ""
if len(sys.argv) == 4:
    region = sys.argv[3]
    region_match = re.fullmatch(r"(\w+):(\d+)-(\d+)", region)
    if region_match is None:
        raise ValueError("Region spec is chr:start-stop")
    seg = region_match.group(1)
    start = int(region_match.group(2))
    stop = int(region_match.group(3))


num_skipped = 0
with open(output_fname, "w") as fout:
    fout.write(
        "\t".join(["site_idx", "alleles", "num_ambigs", "ambig_alleles", "lvl1_site"])
        + "\n"
    )
    for idx, site in enumerate(jvcf["Sites"]):
        if region != "":
            pos = site["POS"]
            if not start <= pos <= stop or seg != site["SEG"]:
                num_skipped += 1
                continue

        num_ambigs = 0
        gtyped_alleles = set()
        alleles = site["ALS"]
        for samp_idx, ft in enumerate(site["FT"]):
            if "AMBIG" in ft:
                num_ambigs += 1
                gt = site["GT"][samp_idx][0]
                if gt is not None:
                    gtyped_alleles.add(alleles[gt])

        num_ambigs = sum(["AMBIG" in ft for ft in site["FT"]])
        lvl1 = "1" if idx in lvl1_sites else "0"
        if len(gtyped_alleles) == 0:
            gtyped_alleles.add(" ")
        fout.write(
            "\t".join(
                [
                    str(idx),
                    ",".join(alleles),
                    str(num_ambigs),
                    ",".join(gtyped_alleles),
                    lvl1,
                ]
            )
            + "\n"
        )

if region != "":
    print(f"Skipped {num_skipped} out of {idx + 1} sites not matching {region}")
