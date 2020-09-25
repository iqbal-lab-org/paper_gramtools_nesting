library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot tb_bigdel alignment stats across tools")
p <- add_argument(p, "input_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "gramtools_commit", help="")

plot_ecdf <- function(dataset, output_dir, commit, name){
  ecdf_plot <- ggplot(dataset, aes(x = NM, colour = condition)) + stat_ecdf(geom="step") + 
    ylab("fraction of sequences") + xlab("edit distance") + 
    theme(text=element_text(size = 13))
  fname <- sprintf("ecdf_plot_%s_%s.pdf", commit, name)
  ggsave(file.path(output_dir,fname),width = 9, height = 6, plot=ecdf_plot)
}

write_mean_eddist <- function(dataset, output_dir, commit, name){
# Write mean edit distance of aligned sequences by condition
  mean_eddist <- dataset %>% group_by(condition) %>% summarise(dist=mean(NM,na.rm=TRUE))
  fname <- sprintf("mean_eddist_%s_%s.txt", commit, name)
  write_tsv(mean_eddist, file.path(output_dir, fname))
}

plot_upset <- function(dataset, output_dir, commit, name){
  # Convert to long format: matrix of 1s and 0s saying whether each seq was found
  df_long_numeric <- dataset %>%
    mutate(found = as.numeric(!is.na(NM))) %>%
    select(sample, gene, condition, found) %>%
    spread(condition, found) %>%
    replace(is.na(.), 0)
  fname <- sprintf("upset_plot_%s_%s.pdf", commit, name)
  pdf(file=file.path(output_dir,fname),width = 9, height = 6, onefile=FALSE)
  plot<-upset(as.data.frame(df_long_numeric), sets=,order.by="freq", text.scale=1.4)
  print(plot)
  dev.off()
}

write_num_aligned <- function(dataset, output_dir, commit, name){
# Write number of aligned sequences by condition
  num_aligned <- dataset %>% filter(!is.na(NM)) %>% 
    group_by(condition) %>% summarise(n())
  fname <- sprintf("num_aligned_%s_%s.txt", commit, name)
  write_tsv(num_aligned, file.path(output_dir, fname))
}


argv <- parse_args(p)
gram_commit = argv$gramtools_commit
outdir = argv$output_dir
dir.create(outdir)

#df_unfiltered <- as_tibble(read.csv("/tmp/stats.tsv", sep="\t"))
df_unfiltered <- as_tibble(read.csv(argv$input_tsv, sep="\t"))

# Remove commit from gramtools condition name
conditions = as.character(unique(df_unfiltered$condition))
gram_cond = conditions[grep("gramtools_*", conditions)]
df_unfiltered$condition = replace(as.vector(df_unfiltered$condition), df_unfiltered$condition == gram_cond, "gramtools")
conditions = as.character(unique(df_unfiltered$condition))

# Filter out alignments with worse MAPQ than baseline_ref seqs and using masks of places where ilmn
# reads map poorly to truth assemblies
df_mapq <- df_unfiltered %>% filter(is.na(MAPQ) | MAPQ > 40)
df_mask <- df_unfiltered %>% filter(mask_overlap <= 0.10)

## Plot empirical cumulative distribution functions ##
plot_ecdf(df_unfiltered, outdir, gram_commit, "unfiltered")
plot_ecdf(df_mapq, outdir, gram_commit, "filtered_mapq40")
plot_ecdf(df_mask, outdir, gram_commit, "filtered_mask10percent")

write_mean_eddist(df_unfiltered, outdir, gram_commit, "unfiltered")
write_mean_eddist(df_mapq, outdir, gram_commit, "filtered_mapq40")


## Plot set intersections: upset plot ##

# Using upsetR package, gives set size as well
library(UpSetR)
# The filter on MAPQ achieves following:
#    - Keeps mapped sequences only (MAPQ is not na)
#    - Keeps mapped sequences where delta_MAPQ is NA: this means we recover a seq that baseline_ref did not
#    - Keeps mapped sequences only if delta_MAPQ is not negative: removes worse alignments

plot_upset(df_unfiltered, outdir, gram_commit, "unfiltered")
plot_upset(df_mapq, outdir, gram_commit, "filtered_mapq40")
plot_upset(df_mask, outdir, gram_commit, "filtered_mask10percent")

write_num_aligned(df_unfiltered, outdir, gram_commit, "unfiltered")
write_num_aligned(df_mapq, outdir, gram_commit, "filtered_mapq40")

## Using ggupset package, interfaces with ggplot
#library(ggupset)
#df_long <- df %>% mutate(found = !is.na(NM)) %>% select(sample, gene, condition, found) %>% spread(condition, found)
#df_upset <- df_long %>% mutate(elem=paste(sample,gene,sep="_")) %>% select(-sample, -gene) %>%
#  gather(Condition, Member, -elem)  %>% filter(Member) %>% select(-Member)
#df_upset <- df_upset %>% group_by(elem) %>% summarise(Sequences=list(Condition))
#
#df_upset %>%  ggplot(aes(x = Sequences)) + geom_bar() + scale_x_upset() + 
#  geom_text(stat='count', aes(label=..count..), vjust = -1)

### Utilities ###
# Show sequences missed by gramtools but found by others
#df_long %>% filter(baseline_ref == TRUE & gramtools == FALSE & vg == TRUE)