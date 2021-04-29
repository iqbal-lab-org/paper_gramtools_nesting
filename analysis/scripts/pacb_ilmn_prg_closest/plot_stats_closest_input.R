library(ggplot2)
library(readr)
library(dplyr)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot performance difference between gramtools and closest input in panel")
p <- add_argument(p, "validation_tsv", help="")
p <- add_argument(p, "closest_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "gramtools_commit", help="")

argv <- parse_args(p)
output_dir <- normalizePath(argv$output_dir)
dir.create(output_dir, recursive=TRUE)

#setwd("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work/jupyter")
#gmtools_commit <- "124321a0"
#closest_df <- readr::read_tsv("closest_stats.tsv")
#validation_df <- readr::read_tsv("validation_stats.tsv")
output_dir <- "."

gmtools_commit <- argv$gramtools_commit

closest_df <- readr::read_tsv(argv$closest_tsv)
closest_df <- rename(closest_df,closest_input_NM=NM)

validation_df <- readr::read_tsv(argv$validation_tsv)

validation_df <- validation_df %>% filter(grepl("gramtools",condition) & grepl(gmtools_commit, condition))
validation_df$condition <- sub("gramtools_.+","gramtools",validation_df$condition)

validation_df <- rename(validation_df,gramtools_NM=NM)
merged_df <- validation_df %>% inner_join(closest_df, by=c("sample","gene"))

p<-ggplot(filter(merged_df, grepl("mapq_20",condition.y)),aes(x=closest_input_NM,y=gramtools_NM)) + geom_count() + 
  scale_size_area("sample count",breaks=c(1,5)) +
  geom_line(aes(x=closest_input_NM,y=closest_input_NM),linetype=2,alpha=0.8) +
  facet_wrap(vars(gene),scales="free") +
  scale_x_continuous("Distance to truth (closest input)") +
  scale_y_continuous("Distance to truth (gramtools)")
ggsave(paste0(output_dir,"/closest_input_mapq_20.pdf"), p,height=10,width=12)

