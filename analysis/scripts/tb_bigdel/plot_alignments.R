library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot tb_bigdel alignment stats across tools")
p <- add_argument(p, "input_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "gramtools_commit", help="")


argv <- parse_args(p)
dir.create(argv$output_dir)

#df <- as_tibble(read.csv("/tmp/stats.tsv", sep="\t"))
df <- as_tibble(read.csv(argv$input_tsv, sep="\t"))

conditions = as.character(unique(df$condition))
gram_cond = conditions[grep("gramtools_*", conditions)]
df$condition = replace(as.vector(df$condition), df$condition == gram_cond, "gramtools")
conditions = as.character(unique(df$condition))

## Plot empirical cumulative distribution function ##
df_filt = df %>% filter(NM < 0.005)
ecdf_plot <- ggplot(df_filt, aes(x = NM, colour = condition)) + stat_ecdf(geom="step") + 
  ylab("fraction of sequences") + xlab("edit distance")
fname <- sprintf("ecdf_plot_%s.pdf", argv$gramtools_commit)
ggsave(file.path(argv$output_dir,fname),width = 9, height = 6, plot=ecdf_plot)

## Plot set intersections: upset plot ##
# Convert to long format


# Using upsetR package, gives set size as well
library(UpSetR)
df_long_numeric <- df %>% mutate(found = as.numeric(!is.na(NM))) %>% select(sample, gene, condition, found) %>% spread(condition, found)
#pdf(file=file.path("/tmp","upset_plot.pdf"),width = 9, height = 6)
fname <- sprintf("upset_plot_%s.pdf", argv$gramtools_commit)
pdf(file=file.path(argv$output_dir,fname),width = 9, height = 6, onefile=FALSE)
upset(as.data.frame(df_long_numeric), sets=,order.by="freq")
dev.off()


## Using ggupset package, interfaces with ggplot
#library(ggupset)
#df_long <- df %>% mutate(found = !is.na(NM)) %>% select(sample, gene, condition, found) %>% spread(condition, found)
#df_upset <- df_long %>% mutate(elem=paste(sample,gene,sep="_")) %>% select(-sample, -gene) %>%
#  gather(Condition, Member, -elem)  %>% filter(Member) %>% select(-Member)
#df_upset <- df_upset %>% group_by(elem) %>% summarise(Sequences=list(Condition))
#
#df_upset %>%  ggplot(aes(x = Sequences)) + geom_bar() + scale_x_upset() + 
#  geom_text(stat='count', aes(label=..count..), vjust = -1)

# Show sequences missed by gramtools but found by others
#df_long %>% filter(baseline_ref == TRUE & gramtools == FALSE & vg == TRUE)
