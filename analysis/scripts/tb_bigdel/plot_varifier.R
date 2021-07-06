library(readr)
library(dplyr)
library(ggplot2)
library(argparser, quietly=TRUE)
library(ggsci)
#library(cowplot)

rm(list=ls())

p <- arg_parser("Plot tb_bigdel varifier-based assessment")
p <- add_argument(p, "varifier_stats_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "deletion_bed", help="Bed file specifying the studied deletion region coordinates")

argv <- parse_args(p)
outdir = argv$output_dir
del_bed <- argv$deletion_bed
dir.create(outdir)
#outdir <- "."
#setwd("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work/varif_perf")
#df <- readr::read_tsv("stats_all.tsv")
#del_bed <- "vars_300.bed"

theme_set(theme(text=element_text(size = 15)))

df <- readr::read_tsv(argv$varifier_stats_tsv)
df <- df %>% mutate(Tool=sub("gramtools_.*","gramtools",Tool))
df <- df %>% mutate(Tool=sub("input_variants","input variants",Tool))
df$Tool <- factor(df$Tool, levels=c("gramtools","graphtyper2","vg","input variants"))

var_type_distribs <- ggplot(df,aes(x=Event_size)) + 
  geom_histogram(aes(fill=Tool),position="dodge") + 
  facet_wrap(vars(Var_type))

write_global_perf <- function(dataset, output_dir, name){
# Write mean recall/precision, across all variant types, for each tool
  global_perfs <- dataset %>% group_by(Tool, Metric) %>% summarise(performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum))
  fname <- sprintf("%s.txt", name)
  write_tsv(global_perfs, file.path(output_dir, fname))
  return(global_perfs)
}

## Classify variants into size bins ##
bin_vartype <- function(var_type, event_size){
 if (var_type %in% c("SNP","MNP")) return(var_type)
  size_name <- ""
 if (event_size <= 10) size_name <- " [1,10]"
  else if (event_size > 50) size_name <- " [51, inf)"
  else size_name <- " [11, 50]"
  return(paste0(var_type,size_name))
}
events <- c("[1,10]","[11, 50]","[51, inf)")
all_bins <- c(paste(rep("DEL"),events), paste(rep("INS"),events),c("SNP","MNP")) 

df <- df %>% rowwise() %>% mutate(binned_vartype=bin_vartype(Var_type, Event_size))

df %>% group_by(Tool, Metric, binned_vartype) %>% summarise(performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum))

## Plot performance of each tool broken down by variant bin
df_perfs <- df %>% group_by(binned_vartype, Tool, Metric) %>% 
  summarise(num_vars=n(), performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum)) %>%
  filter(Tool != "input variants")

by_var_and_tool <- ggplot(df_perfs,aes(x=Metric,y=performance, fill=Tool)) + 
  geom_col(position="dodge") + 
  scale_y_continuous("Performance") +
  facet_wrap(vars(binned_vartype)) + 
  geom_text(aes(label=num_vars), position = position_dodge(width = 1),size=5,vjust=-0.05) +
  theme(text=element_text(size=17))
ggsave(file.path(outdir,"breakdown_by_var_and_tool.pdf"),by_var_and_tool,height=11,width=14)

## The recall is suprisingly low, especially for SNPs. 
global_perfs <- write_global_perf(df, outdir,"mean_perfs_all_vars")
global_perfs

## Plot recall/precision as above but restricting to input variants for the analysed samples which
## are not missing in the truth assemblies; ie, variants expected to be called by genome graph callers.

# Confirm that group by captures 4 tools per group
#test <- df %>% filter(Metric == "recall") %>% group_by(POS, Var_type, Sample) %>% count()
#table(test$n)

no_input_missing <- df %>% group_by(POS, Var_type, Sample, Metric) %>% 
  filter(! any(Eddist_perf_num!=Eddist_perf_denum & Tool == "input variants" & Metric == "recall"))
df_perfs <- no_input_missing %>% group_by(binned_vartype, Tool, Metric) %>% 
  summarise(num_vars=n(), performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum)) %>%
  filter(Tool != "input variants")
by_var_and_tool_expected_from_inputs<-ggplot(df_perfs,aes(x=Metric,y=performance, fill=Tool)) + 
  geom_col(position="dodge") + 
  scale_y_continuous("Performance") +
  facet_wrap(vars(binned_vartype)) +
  geom_text(aes(label=num_vars), position = position_dodge(width = 1),size=5,vjust=-0.05)+
  theme(text=element_text(size=17))
ggsave(file.path(outdir,"breakdown_by_var_and_tool_expected_from_VCF_input.pdf"),
       by_var_and_tool_expected_from_inputs,height=11,width=14)


## That's better:
global_perfs_no_input_missing <- write_global_perf(no_input_missing, outdir, "mean_perfs_all_vars_in_input")
global_perfs_no_input_missing

## We thus have evidence that missing recall in the genome graph tools is likely due to the variants being
## missing in the input VCFs, ie missed by the input variant caller(s) [clockwork]

## We can visualise this is another way. Plot missing recall along the genome, to locate hotspots and
## to see if large missingness co-occurs between the genome graph callers and the input variants
missing_call_df <- filter(df, Metric == "recall" & Eddist_perf_num == 0)
# To plot relative frequency: geom_histogram(aes(y=stat(count)/sum(count)))
p<-ggplot(missing_call_df,aes(x=POS,fill=Tool)) + 
  geom_histogram(bins=100) + 
  facet_wrap(vars(Var_type),ncol=1,scales="free")


deletion_regions <- readr::read_tsv(del_bed,col_names=c("contig","start","stop","name")) %>%
  mutate(midpoint = as.integer((start + stop)/2))
p <- p + geom_vline(data=deletion_regions,aes(xintercept=midpoint),linetype="dotted",colour="salmon3",size=0.8) + 
  scale_y_continuous("Missed variant count",trans="log10") + scale_x_continuous("Genomic position") +
  theme(text=element_text(size=17))
ggsave(file.path(outdir,"missing_recall_breakdown.pdf"),p, 
       height=11,width=14)
#p + scale_x_continuous(breaks=deletion_regions$midpoint,minor_breaks=NULL, labels=NULL)

## How much gain/loss through using a genome graph? (which has other sample calls in it than those analysed)
impact_df <- df %>% filter(Metric == "recall") %>%
  mutate(Eddist_perf_num = max(Eddist_perf_num,0)) %>%
  mutate(recall = Eddist_perf_num / Eddist_perf_denum) %>%
  group_by(POS, Var_type, Sample) %>% 
  mutate(genome_graph_impact = recall - recall[which(Tool == "input variants")])

#my_palette = c("#0072B2","#009E73","#F3730F","#E84BA2","#E69F00","#56B4E9")
# to use manual palette: + scale_fill_manual(values=my_palette)
#show_col(pal_lancet()(7))
impact_distrib <- ggplot(filter(impact_df, Tool != "input variants"), aes(x=genome_graph_impact, fill=Tool)) + 
  facet_wrap(vars(Var_type), ncol=1,scale="free") + 
  scale_x_continuous("Fraction recall gain/loss compared to input variants") + 
  scale_y_continuous("Number of variants") + 
  geom_histogram(bins=22,position="dodge",alpha=0.8) + 
  scale_fill_lancet()

impact_means <- impact_df %>% group_by(Var_type,Tool) %>% summarise(mean=mean(genome_graph_impact))
impact_distrib <- impact_distrib + geom_vline(data=filter(impact_means, Tool != "input variants"), 
                            aes(xintercept=mean,colour=Tool), linetype=5, size=0.65, alpha=0.9) +
  scale_colour_lancet()

ggsave(file.path(outdir,"genome_graph_impact.pdf"),impact_distrib, height=10,width=12)
