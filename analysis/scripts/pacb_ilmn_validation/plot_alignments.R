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
width=11
height=7

# turn condition into a factor and reorder it, so that when we facet_wrap the plots,
# we get the 3D7-based plots on the top row, and the personalised ref plots on the 2nd.
df2 <-  mutate(df, condition=fct_relevel(as.factor(condition), "baseline_ref","samtools_baseline_ref", "cortex_baseline_ref",gram_commit, "samtools_pers_ref", "cortex_pers_ref"))
condition.labs=c("3D7", "samtools(3D7)", "cortex(3D7)", "PR", "samtools(PR)", "cortex(PR)")
names(condition.labs)=c("baseline_ref", "samtools_baseline_ref", "cortex_baseline_ref",gram_commit, "samtools_pers_ref", "cortex_pers_ref")

ggplot(data = filter(df2, gene=="EBA175"),  aes(x=NM))+geom_histogram(binwidth=0.005, color = "black", fill = "lightblue", alpha=0.7) + scale_y_continuous(breaks=seq(0,14,1))+ facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +theme_light()+ xlab("scaled edit distance")+ theme(axis.text.x = element_text(angle = 45))+ expand_limits(x=c(0,0.15))
ggsave(file.path(argv$outdir, "R_NM_EBA175.pdf"), width=width, height=height)

ggplot(data = filter(df2, gene=="AMA1"),  aes(x=NM))+geom_histogram(binwidth=0.005, color = "black", fill = "orange", alpha=0.7) + scale_y_continuous(breaks=seq(0,14,1))+ facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +theme_light()+ xlab("scaled edit distance")+ theme(axis.text.x = element_text(angle = 45))+ expand_limits(x=c(0,0.15))
ggsave(file.path(argv$outdir, "R_NM_AMA1.pdf"), width=width, height=height)

ggplot(data = filter(df2, gene=="DBLMSP"),  aes(x=NM))+geom_histogram(binwidth=0.005, color = "black", fill = "green", alpha=0.7) + scale_y_continuous(breaks=seq(0,14,1))+ facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +theme_light()+ xlab("scaled edit distance")+ theme(axis.text.x = element_text(angle = 45))+ expand_limits(x=c(0,0.15))
ggsave(file.path(argv$outdir, "R_NM_DBLMSP.pdf"), width=width, height=height)

ggplot(data = filter(df2, gene=="DBLMSP2"),  aes(x=NM))+geom_histogram(binwidth=0.005, color = "black", fill = "darkred", alpha=0.7) + scale_y_continuous(breaks=seq(0,14,1))+ facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +theme_light()+ xlab("scaled edit distance")+ theme(axis.text.x = element_text(angle = 45))+ expand_limits(x=c(0,0.15))
ggsave(file.path(argv$outdir, "R_NM_DBLMSP2.pdf"), width=width, height=height)

# For going into main paper, don't expand limits
ggplot(data = filter(df2, gene=="DBLMSP2"),  aes(x=NM))+geom_histogram(binwidth=0.005, color = "black", fill = "darkred", alpha=0.7) + scale_y_continuous(breaks=seq(0,14,1))+ facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +theme_light()+ xlab("scaled edit distance")+ theme(axis.text.x = element_text(angle = 45))
ggsave(file.path(argv$outdir, "R_NM_DBLMSP2_noexpand.pdf"), width=width, height=height)
