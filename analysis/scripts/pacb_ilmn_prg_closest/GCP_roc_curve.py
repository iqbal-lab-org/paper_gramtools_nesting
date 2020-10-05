"""
Computes a precision/recall ROC curve from a tsv describing truth vs genotyped jvcf results.
Only analyses:
    - Regions where mapper could map the sequence to truth assembly, both for gramtools and for closest in prg
    - Regions where this mapping has MAPQ > 40, to be reasonably confident it's placed where it should be in truth assembly
    - Regions (sample x genes) where gramtools has found a NM that is same or worse than closest in prg, otherwise the 'truth' jvcf is not the best gramtools could do (it is worse than it did)
    - Sites which do not have children, to avoid double-counting them
"""
from typing import NewType, Optional, Set, Dict, List
from pathlib import Path
from dataclasses import dataclass
from math import isnan

import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

classif_fields = set(["TP", "TN", "FP", "FN"])
Classif = NewType("Classif", str)


@dataclass
class EvaluatedGene:
    gramtools_NM: float
    prg_closest_NM: Optional[float] = None
    delta_prg_closest: Optional[float] = None

    def delta(self) -> Optional[float]:
        if self.prg_closest_NM is None:
            raise ValueError(
                f"prg_closest_NM attribute is None, it was not set during parsing of prg_closest tsv"
            )
        if self.delta_prg_closest is None:
            self.delta_prg_closest = self.gramtools_NM - self.prg_closest_NM
        return self.delta_prg_closest

    def set_prg_closest(self, value: float) -> None:
        if self.prg_closest_NM is not None:
            raise ValueError(
                f"Cannot set prg_closest_NM: already has value {self.prg_closest_NM}. Duplicate in prg_closest tsv?"
            )
        self.prg_closest_NM = value


Deltas = Dict[str, EvaluatedGene]


def load_gramtools_tsv(gramtools_tsv) -> Deltas:
    deltas = dict()
    df = pd.read_table(gramtools_tsv)
    num_no_NM, num_low_mapq = 0, 0
    for row in df.itertuples():
        if not row.condition.startswith("gramtools"):
            continue
        if isnan(row.NM):
            num_no_NM += 1
            continue
        if row.MAPQ <= 40:
            num_low_mapq += 1
            continue
        key = f"{row.sample}_{row.gene}"
        if key in deltas:
            raise ValueError(f"region {key} found twice in the gramtools tsv")
        deltas[key] = EvaluatedGene(row.NM)
    print(
        f"Gramtools genotyping \n Num no NM: {num_no_NM} \n Num MAPQ <= 40: {num_low_mapq}"
        f"\n Num used regions: {len(deltas)}"
    )
    return deltas


def load_prg_closest_tsv(prg_closest_tsv, deltas: Deltas) -> None:
    num_no_NM = 0
    df = pd.read_table(prg_closest_tsv)
    for row in df.itertuples():
        if not "mapq_40" in row.condition:
            continue
        key = f"{row.sample}_{row.gene}"
        if isnan(row.NM):
            num_no_NM += 1
            if key in deltas:
                # Gramtools could map it, but no input seq in prg is mappable, so remove region from evaluation
                deltas.pop(key)
            continue
        if key not in deltas:
            continue
        deltas[key].set_prg_closest(row.NM)
    print(f"Prg_closest \n Num no NM: {num_no_NM}")


class EvaluatedSite:
    GCP: float
    classif: Classif

    def __init__(self, tsv_line: pd.Series):
        self.GCP = tsv_line.GCP
        classif = tsv_line.classif
        if classif not in classif_fields:
            raise ValueError(f"{classif} not in {classif_fields}")
        self.classif = classif


Classifs = List[EvaluatedSite]


def load_eval_tsv(eval_tsv, deltas: Deltas) -> Classifs:
    classifs = list()
    df = pd.read_table(eval_tsv)
    for row in df.itertuples():
        key = f"{row.sample}_{row.gene}"
        if key not in deltas:
            continue
        if deltas[key].delta() < 0:
            continue
        if row.num_child_sites != 0:
            continue
        if row.ambiguous:
            continue
        classifs.append(EvaluatedSite(row))
    classifs.sort(key=lambda elem: elem.GCP, reverse=True)
    return classifs


ClassifCounts = Dict[Classif, int]


def get_fpr(counts: ClassifCounts) -> float:
    TNs = counts["TN"]
    FPs = counts["FP"]
    return FPs / (TNs + FPs)


def get_precision(counts: ClassifCounts) -> float:
    TPs = counts["TP"]
    FPs = counts["FP"]
    if (TPs + FPs) == 0:
        return 0
    return TPs / (TPs + FPs)


def get_tpr(counts: ClassifCounts) -> float:
    TPs = counts["TP"]
    FNs = counts["FN"]
    return TPs / (TPs + FNs)


def count_classif(sites: Classifs, target_classifs: Set[Classif]) -> int:
    assert classif_fields.issuperset(target_classifs)
    result = 0
    for site in sites:
        if site.classif in target_classifs:
            result += 1
    return result


def check_no_change(call_rates: Dict) -> bool:
    same_last_two = lambda list_elem: list_elem[-1] == list_elem[-2]
    return all(map(same_last_two, call_rates.values()))


def get_roc_values(sites: Classifs) -> pd.DataFrame:
    NUM_STEPS = 30
    tprs, fprs = list(), list()
    counts = {elem: 0 for elem in classif_fields}
    # These two counts will be constant denominators to TPR and FPR respectively
    num_true_calls = count_classif(sites, {"TP", "FN"})
    num_no_calls = count_classif(sites, {"TN", "FP"})
    # First data point: pretend like all sites have been filtered out and called null.
    counts["FN"] = num_true_calls
    counts["TN"] = num_no_calls
    call_rates = {
        "tpr": [get_tpr(counts)],
        "fpr": [get_fpr(counts)],
        "one_minus_precision": [get_precision(counts)],
    }

    step = 0
    step_size = (
        len(sites) // NUM_STEPS
    )  # steps equally spaced in the distribution of GCP
    for _ in range(NUM_STEPS):
        subarray = sites[step : step + step_size]
        TPs_in_range = count_classif(subarray, {"TP"})
        counts["TP"] += TPs_in_range
        counts["FN"] -= TPs_in_range

        FPs_in_range = count_classif(subarray, {"FP"})
        counts["FP"] += FPs_in_range
        counts["TN"] -= FPs_in_range
        call_rates["tpr"].append(get_tpr(counts))
        call_rates["fpr"].append(get_fpr(counts))
        call_rates["one_minus_precision"].append(1 - get_precision(counts))
        step += step_size

        if check_no_change(call_rates):
            break

    result = pd.DataFrame(call_rates)
    return result


@click.command()
@click.argument("gramtools_tsv", type=click.Path(exists=True))
@click.argument("closest_tsv", type=click.Path(exists=True))
@click.argument("evaluation_tsv", type=click.Path(exists=True))
@click.argument("output_dir", type=Path)
def main(gramtools_tsv, closest_tsv, evaluation_tsv, output_dir):
    output_dir.mkdir(exist_ok=True)

    deltas: Dict[str, EvaluatedGene] = load_gramtools_tsv(gramtools_tsv)
    load_prg_closest_tsv(closest_tsv, deltas)
    num_nonneg_delta = 0
    # Compute the deltas and fail early if there is a region with NM in gramtools but not in prg_closest
    for key, elem in deltas.items():
        try:
            delta = elem.delta()
            if delta >= 0:
                num_nonneg_delta += 1
        except ValueError as err:
            raise err
    print(
        f"Total num loaded regions with NM in both gramtools and prg_closest: {len(deltas)}"
    )
    print(
        f"Total num loaded regions with NM in gramtools >= NM in prg_closest: {num_nonneg_delta}"
    )
    classifs = load_eval_tsv(evaluation_tsv, deltas)
    print(f"Num evaluated sites: {len(classifs)}")
    roc_df = get_roc_values(classifs)
    roc_df.to_csv(output_dir / "ROC_stats.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
