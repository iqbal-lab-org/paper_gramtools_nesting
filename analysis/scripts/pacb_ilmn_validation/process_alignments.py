import sys
from collections import namedtuple
from typing import Dict, Tuple
from pathlib import Path

from pysam import AlignmentFile
import seaborn as sns
import pandas as pd

def usage():
    print(f'usage: {sys.argv[0]} input_dir output_dir [input_bed]\n'
            'The input_dir should contain the .sam files to analyse')
    exit(1)

class MultipleAlignmentsError(Exception):
    pass

Scores = namedtuple('Scores', ['NM', 'AS'])
GeneScores = Dict[str, Scores]
BaselineScores: Dict[str, GeneScores] = dict()

def get_sample_and_condition_name(sam_fname: Path) -> Tuple[str, str]:
    elements = sam_fname.stem.split("_")
    sample = elements[-1]
    condition = "_".join(elements[:-1])
    return sample, condition

def load_gene_lengths(bed_fname: Path) -> Dict[str, int]:
    result = dict()
    with bed_fname.open() as fin:
        for line in fin:
            cols = line.split("\t")
            result[cols[3].strip()] = int(cols[2]) - int(cols[1]) + 1
    return result

def get_scores(sam_fname: Path, gene_lengths: Dict[str,int]) -> GeneScores:
    result: GeneScores = dict()
    samfile = AlignmentFile(sam_fname, "r")
    for read in samfile.fetch(until_eof=True):
        gene_name = read.query_name
        if gene_name in result:
            raise MultipleAlignmentsError(f"Several alignments to gene {gene_name} in file {sam_fname}")
        NM = read.get_tag('NM')
        if len(gene_lengths) > 0:
            NM = NM / gene_lengths[gene_name]
        AS = read.get_tag('AS') * -1 # Convert to positive so that high is bad, low is good, like for NM.
        result[gene_name] = Scores(NM, AS)
    return result

def compute_delta(scores_1: Scores, scores_2: Scores) -> Scores:
    newattrs = dict()
    for attr in Scores._fields:
        newattrs[attr] = getattr(scores_2, attr) - getattr(scores_1, attr)
    return Scores(**newattrs)

def get_delta_scores(baseline_genes: GeneScores, query_genes: GeneScores) -> GeneScores:
    """Computes the scores of genes in `query_genes` minus the scores of genes in `baseline_genes`"""
    if baseline_genes.keys() != query_genes.keys():
        raise ValueError(f"Cannot get score deltas for dicts without same keys: {baseline_genes} vs {query_genes}")
    result: GeneScores = dict()
    for gene_name, scores_1 in baseline_genes.items():
        scores_2 = query_genes[gene_name]
        result[gene_name] = compute_delta(scores_1, scores_2)
    return result

if __name__ == "__main__":
    if not 3 <= len(sys.argv) <= 4:
        usage()

    input_dir = Path(sys.argv[1]).resolve()
    if not input_dir.exists():
        print(f'Error: {input_dir} not found')
        usage()

    output_dir = Path(sys.argv[2]).resolve()
    output_dir.mkdir(exist_ok=True)
    output_stats = output_dir / "stats.tsv"

    gene_lengths = dict()
    if len(sys.argv) == 4:
        input_bed = Path(sys.argv[3]).resolve()
        if not input_bed.exists():
            print(f'Error: {input_bed} not found')
            usage()
        gene_lengths = load_gene_lengths(input_bed)

    sam_file_list = list(input_dir.glob(f'*.sam'))
    if len(sam_file_list) == 0:
        print(f'Error: no .sam files in {input_dir}')
        usage()

    ### Write stats ####

    # First pass: get baseline scores
    baseline_scores: BaselineScores = dict()
    for sam_fname in sam_file_list:
        sample, condition = get_sample_and_condition_name(sam_fname)
        if condition == "3d7":
           baseline_scores[sample] = get_scores(sam_fname, gene_lengths) 

    # Second pass: get scores and deltas relative to baseline
    with output_stats.open("w") as stats_file:
        fieldnames = ["sample", "gene", "condition", "NM", "AS", "delta_NM", "delta_AS"]
        stats_file.write("\t".join(fieldnames) + "\n")

        for sam_fname in sam_file_list:
            sample, condition = get_sample_and_condition_name(sam_fname)
            scores = get_scores(sam_fname, gene_lengths)
            delta_scores = get_delta_scores(baseline_scores[sample], scores)
            for gene in delta_scores:
                row = [sample, gene, condition, scores[gene].NM, scores[gene].AS, delta_scores[gene].NM, delta_scores[gene].AS]
                stats_file.write("\t".join(map(str,row)) + "\n")


    ### Write plots ###
    stats_data = pd.read_table(str(output_stats), sep="\t")

    NM = sns.catplot(x='condition', y='NM', kind='box', col='gene',data=stats_data)
    NM.savefig(str(output_dir / "NM.pdf"))

    delta_NM = sns.catplot(x='condition', y='delta_NM', kind='box', col='gene',data=stats_data)
    delta_NM.savefig(str(output_dir / "delta_NM.pdf"))

    delta_AS = sns.catplot(x='condition', y='delta_AS', kind='box', col='gene',data=stats_data)
    delta_AS.savefig(str(output_dir / "delta_AS.pdf"))


    ## Analyse stats ##
    ## Get the distribution of NM change between gramtools genotype and cortex_pers_ref
    data = pd.read_csv(str(output_stats), sep="\t")
    def get_NMs(condition): 
      idx1, idx2 = 0, 0 
      for els in condition.iteritems(): 
         if els[1] == 'cortex_pers_ref': 
             idx1 = els[0] 
         if els[1] == 'gramtools_genotype': 
             idx2 = els[0] 
      return data.loc[idx1]["NM"] - data.loc[idx2]["NM"] 
    NM_changes = data.groupby(["sample","gene"])["condition"].agg(get_NMs) 
    NM_changes.to_csv(output_dir / "discov_genotype_NM_diffs.tsv", sep="\t")

