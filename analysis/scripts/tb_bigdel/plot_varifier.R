library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work/varif_perf")


df <- readr::read_tsv("stats_all.tsv")
df <- df %>% mutate(Tool=sub("gramtools_.*","gramtools",Tool))
df$Tool <- factor(df$Tool, levels=c("gramtools","graphtyper2","vg","input_variants"))

var_type_distribs <- ggplot(df,aes(x=Event_size)) + 
  geom_histogram(aes(fill=Tool),position="dodge") + 
  facet_wrap(vars(Var_type))


## Classify variants into size bins ##
bin_vartype <- function(var_type, event_size){
 if (var_type %in% c("SNP","MNP")) return(var_type)
  size_name <- ""
 if (event_size <= 10) size_name <- " [0,10]"
  else if (event_size > 50) size_name <- " [51, inf]"
  else size_name <- " [11, 50]"
  return(paste0(var_type,size_name))
}

all_bins <- c("DEL [0,10]","DEL [11, 50]", "DEL [51, inf]", 
               "INS [0,10]","INS [11, 50]", "INS [51, inf]",
               "SNP", "MNP")

df <- df %>% rowwise() %>% mutate(binned_vartype=bin_vartype(Var_type, Event_size))
df$binned_vartype <- factor(df$binned_vartype, levels=all_bins)


## Plot performance of each tool broken down by variant bin
df_perfs <- df %>% group_by(binned_vartype, Tool, Metric) %>% 
  summarise(num_vars=n(), performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum)) %>%
  filter(Tool != "input_variants")

by_var_and_tool <- ggplot(df_perfs,aes(x=Metric,y=performance, fill=Tool)) + 
  geom_col(position="dodge") + 
  geom_text(aes(label=num_vars), position = position_dodge(width = 1),vjust=-0.25) +
  facet_wrap(vars(binned_vartype))
ggsave("breakdown_by_var_and_tool.pdf",by_var_and_tool,height=10,width=13)

## The recall is suprisingly low, especially for SNPs. 
global_perfs <- df %>% group_by(Tool, Metric) %>% summarise(performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum))
global_perfs

## Plot recall/precision as above but restricting to input variants for the analysed samples which
## are not missing in the truth assemblies; ie, variants expected to be called by genome graph callers.
no_input_missing <- df %>% group_by(POS, Var_type) %>% 
  filter(! any(Eddist_perf_num==0 & Tool == "input_variants" & Metric == "recall"))
df_perfs <- no_input_missing %>% group_by(binned_vartype, Tool, Metric) %>% 
  summarise(num_vars=n(), performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum))
by_var_and_tool_expected_from_inputs<-ggplot(df_perfs,aes(x=Metric,y=performance, fill=Tool)) + 
  geom_col(position="dodge") + 
  geom_text(aes(label=num_vars), position = position_dodge(width = 1),vjust=-0.25) +
  facet_wrap(vars(binned_vartype))
ggsave("breakdown_by_var_and_tool_expected_from_VCF_input.pdf",by_var_and_tool_expected_from_inputs,height=10,width=13)

## That's better:
no_input_missing %>% group_by(Tool, Metric) %>% summarise(performance=sum(Eddist_perf_num) / sum(Eddist_perf_denum))

## We thus have evidence that missing recall in the genome graph tools is likely due to the variants being
## missing in the input VCFs, ie missed by the input variant caller(s) [clockwork]

## We can visualise this is another way. Plot missing recall along the genome, to locate hotspots and
## to see if large missingness co-occurs between the genome graph callers and the input variants
missing_call_df <- filter(df, Metric == "recall" & Eddist_perf_num == 0)
# To plot relative frequency: geom_histogram(aes(y=stat(count)/sum(count)))
p<-ggplot(missing_call_df,aes(x=POS,fill=Tool)) + 
  geom_histogram(bins=300) + 
  facet_wrap(vars(Var_type),ncol=1,scales="free")

deletion_regions <- readr::read_tsv("vars_300.bed",col_names=c("contig","start","stop","name")) %>%
  mutate(midpoint = as.integer((start + stop)/2))
p <- p + geom_vline(data=deletion_regions,aes(xintercept=midpoint),linetype="dotted",colour="salmon3") + 
  scale_y_continuous("Missed variant count",trans="log10") + scale_x_continuous("Genomic position")
ggsave("missing_recall_breakdown.pdf",p, height=10,width=13)
#p + scale_x_continuous(breaks=deletion_regions$midpoint,minor_breaks=NULL, labels=NULL)

test_dels <- df %>% filter(binned_vartype == "DEL>10" & Metric == "recall")
test_snps <- df %>% filter(binned_vartype == "SNP" & Metric == "recall")


