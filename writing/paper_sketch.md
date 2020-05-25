# Paper Contents

## Key message

* Genome graphs with nested variation occur and are biologically relevant
* Mechanisms are needed to perform genotyping on nested genome graphs
* We show a model and algorithm to achieve this
* We show how it can be useful in practice


## Methods

* Gramtools genotyping model
	* Coverage: per base and per-allele with equivalence classes
	* Sequencing errors
* Applying the model to bubbles in the graph recursively
	* Use a total ordering of bubbles based on REF-based POS in genome graph to genotype nested bubbles before their parents
	* Parent bubble considers only calls from its children, not all possible calls
	* Invalidation process to maintain coherence: children of haplogroups called as absent must be called as absent too
* A useful output
	* Hard to see how VCF can express nested variation. Would require multiple references/ overlapping records
	* Made a JSON format that produces overlapping records and a map of the nested variation.


## Results

* Analyse allelic dimorphism/balancing selection in one Pf gene
* Analyse all deletions with variants under them in TB

# Work required

## Genotyping model improvements and evaluation

Need a framework for evaluating genotyping model quality so as to test effect of things like:
* prob distribution used in likelihood terms
* genotype confidence percentiles 
* using different haplogroups only for genotype confidence
* having a genotype that is 'none of the alleles at this site'

Perhaps construct PRGs of regions in TB and Pf, simulate paths (gramtools `simulate`), reads, genotype and compare to truth
from `simulate`.

## Code to process the JSON output

To prove its a rational and useful output and that VCF alone would not be enough.

## Analyses

* Show we get good calls using truth assemblies
* Show the calls are better/more useful than in a non-nested VCF (which is also produced by `gramtools genotype`)
