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

#df = read_tsv("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work/124321a0/stats.tsv")
#gram_commit <- "gramtools_124321a0"
#outdir <- "."
outdir <- argv$outdir

theme_set(theme_classic(base_size = 20))
# turn condition into a factor and reorder it, so that when we facet_wrap the plots,
# we get the 3D7-based plots on the top row, and the personalised ref plots on the 2nd.
df2 <-  mutate(df, condition=fct_relevel(as.factor(condition), "baseline_ref","samtools_baseline_ref", "cortex_baseline_ref",gram_commit, "samtools_pers_ref", "cortex_pers_ref"))
condition.labs=c("3D7", "SAMtools(3D7)", "Cortex(3D7)", "PR", "SAMtools(PR)", "Cortex(PR)")
names(condition.labs)=c("baseline_ref", "samtools_baseline_ref", "cortex_baseline_ref",gram_commit, "samtools_pers_ref", "cortex_pers_ref")

for (gene_name in unique(df2$gene)){
filtered_df <- filter(df2, gene==gene_name)
means <- filtered_df %>% group_by(condition) %>% summarise(mean=mean(NM))
ggplot(filtered_df,  aes(x=NM))+
  geom_histogram(binwidth=0.005, color = "black", fill = "salmon3", alpha=0.7) + 
  scale_y_continuous(breaks=seq(0,14,1))+ 
  facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +
  theme_light()+ 
  xlab("scaled edit distance")+ 
  expand_limits(x=c(0,0.15)) + 
  theme(text = element_text(size=17), axis.text.x = element_text(angle=20)) +
  geom_vline(data=means,aes(xintercept=mean),linetype="dotted",size=0.8)
ggsave(file.path(outdir, paste0("/R_NM_",gene_name,".pdf")), width=width, height=height)
}
p<-ggplot(filtered_df,  aes(x=NM))+
  geom_histogram(bins=20, color = "black", fill = "salmon3", alpha=0.7) + 
  scale_y_continuous(breaks=seq(0,14,1))+ 
  facet_wrap(~condition, labeller=labeller(condition=condition.labs)) +
  theme_light()+ 
  xlab("scaled edit distance")+ 
  #expand_limits(x=c(0,0.15)) + 
  theme(text = element_text(size=18), axis.text.x = element_text(angle=20))
p2<-p + geom_vline(data=means,aes(xintercept=mean),linetype="dotted",size=0.8,show.legend=TRUE)