#!/usr/bin/env python
# coding: utf-8


from typing import List

import pandas as pd
from scipy.cluster.hierarchy import linkage, to_tree, ClusterNode
from scipy.stats import entropy
import json
import click
from collections import Counter, defaultdict
from seaborn import clustermap, color_palette
import seaborn as sns
import matplotlib.pyplot as plt

from jvcf_processing import (
    Region,
    click_get_region,
    is_in_region,
    first_idx_in_region,
)
from common import get_partition, allelic_distinguishability, allelic_specificity


groups=get_partition("heatmaps/DBL_DBLMSP2_hapgs.tsv")


# In[36]:


## Get sample -> country and sample -> dimorphic form hash maps
metadata = pd.read_csv("pf3k_release_5.tsv",sep="\t")
sample_to_country = dict(zip(list(metadata["sample"]),list(metadata["country"])))

sample_to_dimorphic = {sample_name: "form1" for sample_name in groups[0]}
sample_to_dimorphic.update({sample_name: "form2" for sample_name in groups[1]})


# In[37]:


not_none = lambda x: x is not None


# In[38]:



def relative_entropy(gts1, gts2) -> float:
    """What is the KL divergence between the two genotype set prob. distributions?"""
    gt1_distrib, gt2_distrib = get_complete_counts(gts1, gts2)
    if sum(gt1_distrib) == 0 or sum(gt2_distrib) == 0:
        return float("inf")
    return entropy(gt1_distrib,gt2_distrib)



# In[39]:


def get_heterozygosity(values):
    total = sum(values)
    if total == 0:
        raise ValueError(f"Counts sum to 0: {values}")
    het = 1 - sum(map(lambda x: (x/total)**2,values))
    return het

def gt_heterozygosity(gts) -> float:
    c=Counter(list(filter(not_none,gts)))
    return get_heterozygosity(c.values())
    
def cross_heterozygosity(values1, values2) -> float:
    total1 = sum(values1)
    total2 = sum(values2)
    if total1 == 0 or total2 == 0:
        raise ValueError("Zero sum total found")
    prob_same_allele = 0
    for elem1, elem2 in zip(values1, values2):
        prob_same_allele += elem1/total1 * elem2/total2
    return 1 - prob_same_allele


# In[40]:


with open("combined.json") as fin:
    jvcf = json.load(fin)


# In[41]:


region=Region("Pf3D7_10_v3",1432803,1434147)
lvl1_sites = set(jvcf["Lvl1_Sites"])
first_idx = first_idx_in_region(jvcf["Sites"], region)

while first_idx not in lvl1_sites:
    first_idx += 1


# In[42]:


country_dimorphism = defaultdict(lambda: defaultdict(int))
## Figure out how many samples in each country are in each form
for sample in jvcf["Samples"]:
    country = sample_to_country[sample["Name"]]
    try:
        form = sample_to_dimorphic[sample["Name"]]
    except KeyError:
        continue
    country_dimorphism[country][form] += 1
for country in country_dimorphism:
    country_dimorphism[country] = list(country_dimorphism[country].values())


# ## Is the dimorphism present across countries?
# 
# Dataset: 706 samples across three countries
# **Dimorphism** is defined as the basal split of hierarchical clustered dendrogram for which the two groups below are in proportion at least 10%/90% (ie, the smallest group cannot be <10% of all samples below the node)
# This leads to a sample partition of size 329, 358 (total 687).

# In[43]:


country_dimorphism


# In[44]:


cross_het = list()
for counts1 in country_dimorphism.values():
    next_row = list()
    for counts2 in country_dimorphism.values():
        next_row.append(cross_heterozygosity(counts1, counts2))
    cross_het.append(next_row)
mask = [[False * 3],[True, False, False], [True, True, False]]
f, ax = plt.subplots(figsize=(7, 5))
ax = sns.heatmap(cross_het, mask=mask, square=True,vmax=.55,vmin=0, cbar_kws={"label":"interform heterozygosity"})
#ax = sns.heatmap(cross_het, square=True,vmax=.55,vmin=0, cbar_kws={"label":"interform heterozygosity"})
ax.set_xticklabels(country_dimorphism.keys())
ax.set_yticklabels(country_dimorphism.keys())


# In[45]:


cross_het


# ## Analysing the relationship between site-level dimorphism and country-wide genotype diversity
# 
# Goals:
# 
#     - Spot markers of the dimorphism
#     - Illustrate what variation nested sites are picking up
#     - Test whether dimorphic sites correlate with genotype diversity

# In[46]:


dimorphism_calls = defaultdict(list)
country_calls = defaultdict(list)
country_dimorphism_counts = defaultdict(list)
site_is_nested = list()

cur_idx = first_idx
num_sites = 0

while is_in_region(jvcf["Sites"][cur_idx], region):
    dimorphism_site_calls = defaultdict(list)
    country_site_calls = defaultdict(list)
    country_dimo_site_counts = defaultdict(lambda: {"form1":0,"form2":0})
    site = jvcf["Sites"][cur_idx]
    gts = [gt[0] for gt in site["GT"]]
    if len(set(filter(not_none,gts))) < 2:
        cur_idx +=1
        continue
    for sample_idx, gt in enumerate(gts):
        sample_name = jvcf["Samples"][sample_idx]["Name"]
        country = sample_to_country[sample_name]
        country_site_calls[country].append(gt)
        try:
            form = sample_to_dimorphic[sample_name]
            dimorphism_site_calls[form].append(gt)
            if gt is not None:
                country_dimo_site_counts[country][form] += 1
        except KeyError:
            pass

    # Skip a site if it has too many null calls
    skip_site = False
    for country_gts in country_site_calls.values():
        if len(list(filter(not_none, country_gts))) < 5:
            skip_site = True
    if skip_site:
        cur_idx += 1
        continue
    for key,val in country_site_calls.items():
        country_calls[key].append(val)
    for key,val in dimorphism_site_calls.items():
        dimorphism_calls[key].append(val)
    for key, val in country_dimo_site_counts.items():
        country_dimorphism_counts[key].append(list(val.values()))
    site_is_nested.append(False if cur_idx in lvl1_sites else True)
    cur_idx += 1
    num_sites += 1


# In[47]:


num_sites


# In[48]:


## Compute genotype heterozygosity for each site in each country
for country in country_calls:
    heterozygosities = map(gt_heterozygosity,country_calls[country])
    country_calls[country] = list(heterozygosities)
    
## Compute inter-form heterozygosity for each site in each country
for country in country_dimorphism_counts:
    heterozygosities = map(get_heterozygosity,country_dimorphism_counts[country])
    country_dimorphism_counts[country] = list(heterozygosities)


# In[49]:


## Compute measures of per-site dimorphism
dimorphism_sensitivity = []
dimorphism_specificity = []
for site_idx in range(num_sites):
    form1_gts = dimorphism_calls["form1"][site_idx]
    form2_gts = dimorphism_calls["form2"][site_idx]
    dimorphism_sensitivity.append(allelic_distinguishability(form1_gts,form2_gts))
    dimorphism_specificity.append(allelic_specificity(form1_gts,form2_gts))


# In[50]:


plt.hist(dimorphism_sensitivity)


# In[51]:


plt.hist(dimorphism_specificity)


# In[25]:


df = pd.DataFrame.from_dict(country_calls, orient="index")


# In[52]:


p=color_palette("magma",as_cmap=True)
#dimo_sensi_colour=[p(val) for val in dimorphism_sensitivity]
dimo_sensi_colour = ["r" if val > 0.8 else "black" for val in dimorphism_sensitivity]
dimo_speci_colour = list()
for val in dimorphism_specificity:
    if val > 0.3 or val < -0.3:
        dimo_speci_colour.append("r")
    else:
        dimo_speci_colour.append("black")
nestedness=["gray" if val else "navajowhite" for val in site_is_nested]
cl_map=clustermap(df, col_cluster=False,row_cluster=False,
                  col_colors=[dimo_speci_colour,nestedness,dimo_sensi_colour],cmap="viridis",
                  cbar_kws={"label":"heterozygosity"})


# In[199]:


get_ipython().run_line_magic('pinfo', 'cl_map.add_legend')


# In[239]:


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# In[245]:


countries = list(df.index)
fig, ax = plt.subplots()
ax.pcolormesh(df)
ax.set_yticks(np.arange(len(countries)) + 0.5)
ax.set_yticklabels(countries)
im = ax.imshow(df)
fig.colorbar(im)
fig.tight_layout()
plt.show()
#plt.colorbar(label="heterozygosity")


# In[53]:


tmp_dict = country_calls.copy()
tmp_dict["dimorphism_sensitivity"] = dimorphism_sensitivity
wide_df = pd.DataFrame(tmp_dict)
long_df = pd.melt(wide_df, id_vars="dimorphism_sensitivity",value_vars=country_calls.keys(),var_name="Country",value_name="heterozygosity")


# In[54]:


long_df


# In[55]:


g=sns.lmplot(y="heterozygosity",x="dimorphism_sensitivity",data=long_df,hue="Country")
g=g.set_ylabels("genotype heterozygosity")


# In[56]:


tmp_dict = country_dimorphism_counts.copy()
tmp_dict["dimorphism_sensitivity"] = dimorphism_sensitivity
wide_df = pd.DataFrame(tmp_dict)
long_df = pd.melt(wide_df, id_vars="dimorphism_sensitivity",value_vars=country_calls.keys(),var_name="Country",value_name="interform heterozygosity")
sns.lmplot(y="interform heterozygosity",x="dimorphism_sensitivity",data=long_df,hue="Country")

