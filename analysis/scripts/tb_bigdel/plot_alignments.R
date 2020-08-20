library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(ggupset)

df <- as_tibble(read.csv("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work/stats/stats.tsv", sep="\t"))
df_filt = df %>% filter(NM < 0.005)
ggplot(df_filt, aes(x = NM, colour = condition)) + stat_ecdf(geom="step")

# Convert to long format
df_long <- df %>% mutate(found = !is.na(NM)) %>% select(sample, gene, condition, found) %>% spread(condition, found)

df_upset <- df_long %>% mutate(elem=paste(sample,gene,sep="_")) %>% select(-sample, -gene) %>%
  gather(Condition, Member, -elem)  %>% filter(Member) %>% select(-Member)
df_upset <- df_upset %>% group_by(elem) %>% summarise(Sequences=list(Condition))

df_upset %>%  ggplot(aes(x = Sequences)) + geom_bar() + scale_x_upset() + 
  geom_text(stat='count', aes(label=..count..), vjust = -1)

# With upsetR package, gives set size as well
library(UpSetR)
df_long_numeric <- df %>% mutate(found = as.numeric(!is.na(NM))) %>% select(sample, gene, condition, found) %>% spread(condition, found)
upset(as.data.frame(df_long_numeric), sets=as.character(unique(df$condition)),order.by="freq")

# Missed by gramtools but found by others
df_long %>% filter(baseline_ref == TRUE & gramtools_80dfb019 == FALSE & vg == TRUE)
