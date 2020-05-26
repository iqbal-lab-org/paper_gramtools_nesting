Workflows
===========

Simulated data
---------------

Plasmodium DBLMSPs
```````````````````

I re-use the prgs constructed in pf_surfants repository.
I take DBLMSP 1 & 2 prgs built from pf3k vcfs with 5kbp flanks either side of the genes. I build the non-variant bit in between
using gramtools' `encode_prg` taking the sequence in between the genes pulled using `samtools faidx` (remember in the bed, the left side is 
0-based so has to be incremented by 1). I concatenate the three binaries using scripts/starting_prg/concat_prg.py in pf_surfants repo.


