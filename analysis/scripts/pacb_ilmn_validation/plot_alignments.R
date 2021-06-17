library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(forcats)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot alignments")
p <- add_argument(p, "stats_file", help="")
p <- add_argument(p, "outdir", help="")
p <- add_argument(p, "gram_commit", help="")

argv <- parse_args(p)
df = read_tsv(argv$stats_file)
gram_commit = sprintf("gramtools_%s",argv$gram_commit)
width=13
height=10
text_size=18

#setwd("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work/124321a0/")
#df = read_tsv("stats.tsv")
#gram_commit <- "gramtools_124321a0"
#outdir <- "."
outdir <- argv$outdir

theme_set(theme_classic(base_size = 20))

## Plots for each gene
# Rename and order the conditions
condition.labs=c("3D7", "SAMtools(3D7)", "Cortex(3D7)", "PR", "SAMtools(PR)", "Cortex(PR)")
df2 <- df
df2$condition <- replace(df$condition, grepl("gramtools_",df2$condition),"gramtools")
df2$condition <- recode(df2$condition,
  baseline_ref="3D7",samtools_baseline_ref="SAMtools(3D7)",cortex_baseline_ref="Cortex(3D7)",
  gramtools="PR", samtools_pers_ref="SAMtools(PR)",cortex_pers_ref="Cortex(PR)"
)
df2$condition <- factor(df2$condition, levels=condition.labs)

for (gene_name in unique(df2$gene)){
  filtered_df <- filter(df2, gene==gene_name)
  means <- filtered_df %>% group_by(condition) %>% summarise(mean=mean(NM))
  ggplot(filtered_df,  aes(x=NM))+
    geom_histogram(binwidth=0.005, color = "black", fill = "salmon3", alpha=0.7) + 
    scale_y_continuous(breaks=seq(0,14,1))+ 
    facet_wrap(~condition, nrow=2) +
    theme_light()+ 
    xlab("Scaled edit distance") + 
    ylab("Sample count") +
    expand_limits(x=c(0,0.15)) + 
    theme(text = element_text(size=text_size), axis.text.x = element_text(angle=20)) +
    geom_vline(data=means,aes(xintercept=mean),linetype="dotted",size=0.8) +
    geom_text(data=means,aes(x=mean+0.012, y = 10, label=round(mean,3)),size=5)
  ggsave(file.path(outdir, paste0("/R_NM_",gene_name,".pdf")), width=width, height=height)
}

## Plot of all genes
# Rename and order the conditions
condition.labs=c("3D7", "SAMtools", "Cortex", "gramtools")
df3 <- df %>% filter(! grepl("pers_ref",condition))
df3$condition <- replace(df3$condition, grepl("gramtools_",df3$condition),"gramtools")
df3$condition <- recode(df3$condition,
  baseline_ref="3D7",samtools_baseline_ref="SAMtools",cortex_baseline_ref="Cortex"
)
df3$condition <- factor(df3$condition, levels=condition.labs)

means <- df3 %>% group_by(condition) %>% summarise(mean=mean(NM))
p<-ggplot(df3,  aes(x=NM))+
  geom_histogram(bins=20, color = "black", fill = "salmon3", alpha=0.7) + 
  #scale_y_continuous(breaks=seq(0,14,1))+ 
  facet_wrap(~condition,nrow=2) +
  theme_light()+ 
  xlab("scaled edit distance")+ 
  xlab("Scaled edit distance") + 
  ylab("Sequence count") +
  theme(text = element_text(size=text_size), axis.text.x = element_text(angle=20)) +
  geom_vline(data=means,aes(xintercept=mean),linetype="dotted",size=0.8,show.legend=TRUE) +
  geom_text(data=means,aes(x=mean+0.012, y = 40, label=round(mean,3)),size=5)

ggsave(file.path(outdir, paste0("all_genes.pdf")), width=width, height=height)
