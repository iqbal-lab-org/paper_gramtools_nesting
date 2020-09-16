library(ggplot2)
library(tibble)
library(dplyr)

df_unfiltered <- as_tibble(read.csv("/tmp/callsunfiltered_stats.tsv", sep="\t"))
# Remove commit from gramtools condition name
conditions = as.character(unique(df_unfiltered$condition))
gram_cond = conditions[grep("gramtools_*", conditions)]
df_unfiltered$condition = replace(as.vector(df_unfiltered$condition), df_unfiltered$condition == gram_cond, "gramtools")
conditions = as.character(unique(df_unfiltered$condition))


prg_closest <- as_tibble(read.csv("/tmp/prg_stats.tsv", sep="\t"))
df_mapq <- df_unfiltered %>% filter(is.na(MAPQ) | MAPQ > 40)
merged <- bind_rows(df_mapq, prg_closest)
merged_long <- merged %>%
  select(sample, gene, condition, NM) %>%
  spread(condition, NM) 

merged_long <- merged_long %>% mutate(gramtools_diff=gramtools - closest_in_prg_mapq_40)

# In these 32 cases, gramtools precision is 2x worse than vg and 5x worse than gtyper2,
# and the closest input has 97% lower distance than found by gramtools
bad_gtypes <- merged_long %>% filter(gramtools_diff > 0.005)
bad_gtypes %>% summarise_all(mean,na.rm=TRUE)

# Removing 137 cases where gramtools cannot recover the best input, gramtools outperforms other tools
# NB this does not mean the best input is in the graph
good_gtypes <- merged_long %>% filter(gramtools_diff <= 0)
good_gtypes %>% summarise_all(mean,na.rm=TRUE)

ggplot(data=merged_long, aes(x=gramtools_diff)) + geom_histogram()
ggplot(data=merged_long, aes(x=closest_in_prg_mapq_40)) + geom_histogram()
ggplot(data=merged_long, aes(x=closest_in_prg_mapq_40)) + stat_ecdf()
