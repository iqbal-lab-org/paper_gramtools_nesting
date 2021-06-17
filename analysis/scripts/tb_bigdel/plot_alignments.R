library(ggplot2)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(ggsci)
library(argparser, quietly=TRUE)


p <- arg_parser("Plot tb_bigdel alignment stats across tools")
p <- add_argument(p, "input_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "gramtools_commit", help="")


plot_ecdf <- function(dataset, output_dir, commit, name, with_unmapped = FALSE){
  if (with_unmapped) {
    # Replace NA's with a little bit more than each tool's max NM
    temp <- dataset %>% group_by(tool) %>% mutate(local_max_NM=max(NM,na.rm=TRUE) + 0.02)
    maxed <- dplyr::coalesce(dataset$NM, temp$local_max_NM)
    plotted_dataset <- dataset %>% mutate(NM = maxed)
  }
  else {
    plotted_dataset = dataset
  }
  
  ecdf_plot <- ggplot(plotted_dataset, aes(x = NM, colour = tool)) + stat_ecdf(geom="step", 
                                                                                    pad = FALSE) + 
    ylab("Fraction of sequences") + xlab("Edit distance") + labs(colour="Tool") +
    theme(text=element_text(size = 13)) +
    scale_colour_lancet()
  fname <- sprintf("ecdf_plot_%s_%s.pdf", commit, name)
  ggsave(file.path(output_dir,fname),width = 9, height = 6, plot=ecdf_plot)
}

write_mean_eddist <- function(dataset, output_dir, commit, name){
# Write mean edit distance of aligned sequences by tool
  mean_eddist <- dataset %>% group_by(tool) %>% summarise(dist=mean(NM,na.rm=TRUE))
  fname <- sprintf("mean_eddist_%s_%s.txt", commit, name)
  write_tsv(mean_eddist, file.path(output_dir, fname))
}

plot_upset <- function(dataset, output_dir, commit, name, plot_missing = FALSE){
  # Convert to long format: matrix of 1s and 0s saying whether each seq was aligned
  if (plot_missing){
    plotted_dataset <- dataset %>%
      mutate(missing = as.numeric(is.na(NM))) %>%
      select(sample, gene, tool, missing) %>%
      spread(tool, missing)
  }
  else{
  plotted_dataset <- dataset %>%
    mutate(found = as.numeric(!is.na(NM))) %>%
    select(sample, gene, tool, found) %>%
    spread(tool, found)
  }
  
  fname <- sprintf("upset_plot_%s_%s.pdf", commit, name)
  pdf(file=file.path(output_dir,fname),width = 9, height = 6, onefile=FALSE)
  plot<-upset(as.data.frame(plotted_dataset), sets=,order.by="freq", text.scale=1.4)
  print(plot)
  dev.off()
}

write_num_aligned <- function(dataset, output_dir, commit, name){
# Write number of aligned sequences by tool
  num_aligned <- dataset %>% filter(!is.na(NM)) %>% 
    group_by(tool) %>% summarise(n())
  fname <- sprintf("num_aligned_%s_%s.txt", commit, name)
  write_tsv(num_aligned, file.path(output_dir, fname))
}


argv <- parse_args(p)
gram_commit = argv$gramtools_commit
outdir = argv$output_dir
dir.create(outdir)

#setwd("/home/brice/Desktop/main_PhD/analyses/nesting_paper/tmp_work")
#df_unfiltered <- readr::read_tsv("124321a0/minimap2/callsfilterpass_stats.tsv")
#gram_commit <- "124321a0"
#outdir <- "."

df_unfiltered <- as_tibble(read.csv(argv$input_tsv, sep="\t"))
df_unfiltered <- df_unfiltered %>% rename(tool = condition)

# Remove commit from gramtools condition name
tools = as.character(unique(df_unfiltered$tool))
gram_cond = tools[grep("gramtools_*", tools)]
tool_names <- c("gramtools", "graphtyper2", "vg", "baseline(H37Rv)")
df_unfiltered$tool = replace(as.vector(df_unfiltered$tool), df_unfiltered$tool == gram_cond, "gramtools")
df_unfiltered$tool = replace(as.vector(df_unfiltered$tool), df_unfiltered$tool == "baseline_ref", tool_names[4])
# Order the tools so that the colours come out in the same order as in varifier evaluation plot
df_unfiltered$tool <- factor(df_unfiltered$tool, levels=tool_names)

# Convert alignment NM below MAPQ threshold to NA
df_mapq30 <- df_unfiltered
df_mapq30$NM[df_unfiltered$MAPQ < 30] <- NA
df_mapq30_capped_NM <- df_mapq30
df_mapq30_capped_NM$NM[df_mapq30$NM > 0.5 & !is.na(df_mapq30$NM)] <- 0.5
df_mapq40 <- df_unfiltered
df_mapq40$NM[df_unfiltered$MAPQ < 40] <- NA


## Plot empirical cumulative distribution functions ##
plot_ecdf(df_unfiltered, outdir, gram_commit, "unfiltered", with_unmapped = TRUE)
plot_ecdf(df_mapq30, outdir, gram_commit, "filtered_mapq30", with_unmapped = TRUE)
plot_ecdf(df_mapq30_capped_NM, outdir, gram_commit, "filtered_mapq30_maxNM0.5", with_unmapped = TRUE)
plot_ecdf(df_mapq40, outdir, gram_commit, "filtered_mapq40", with_unmapped = TRUE)

write_mean_eddist(df_unfiltered, outdir, gram_commit, "unfiltered")
write_mean_eddist(df_mapq30, outdir, gram_commit, "filtered_mapq30")
write_mean_eddist(df_mapq30_capped_NM, outdir, gram_commit, "filtered_mapq30_maxNM0.5")
write_mean_eddist(df_mapq40, outdir, gram_commit, "filtered_mapq40")


## Plot set intersections: upset plot ##

# Using upsetR package, gives set size as well
library(UpSetR)
# The filter on MAPQ achieves following:
#    - Keeps mapped sequences only (MAPQ is not na)
#    - Keeps mapped sequences where delta_MAPQ is NA: this means we recover a seq that baseline_ref did not
#    - Keeps mapped sequences only if delta_MAPQ is not negative: removes worse alignments

plot_upset(df_unfiltered, outdir, gram_commit, "unfiltered_unmapped", plot_missing = TRUE)
plot_upset(df_mapq30, outdir, gram_commit, "filtered_mapq30_unmapped", plot_missing = TRUE)
plot_upset(df_mapq40, outdir, gram_commit, "filtered_mapq40_unmapped", plot_missing = TRUE)

write_num_aligned(df_unfiltered, outdir, gram_commit, "unfiltered")
write_num_aligned(df_mapq30, outdir, gram_commit, "filtered_mapq30")
write_num_aligned(df_mapq40, outdir, gram_commit, "filtered_mapq40")

## Using ggupset package, interfaces with ggplot
#library(ggupset)
#df_long <- df %>% mutate(found = !is.na(NM)) %>% select(sample, gene, tool, found) %>% spread(tool, found)
#df_upset <- df_long %>% mutate(elem=paste(sample,gene,sep="_")) %>% select(-sample, -gene) %>%
#  gather(Condition, Member, -elem)  %>% filter(Member) %>% select(-Member)
#df_upset <- df_upset %>% group_by(elem) %>% summarise(Sequences=list(Condition))
#
#df_upset %>%  ggplot(aes(x = Sequences)) + geom_bar() + scale_x_upset() + 
#  geom_text(stat='count', aes(label=..count..), vjust = -1)

### Utilities/Debugging ###
# Show sequences missed by gramtools but found by others
#df_long %>% filter(baseline_ref == TRUE & gramtools == FALSE & vg == TRUE)

### Minimap2 results: look at unmapped, low mapq, and gramtools high edit distance sequences ##
#df_unfiltered_minimap2 <- as_tibble(read.csv("/home/brice/Desktop/main_PhD/analyses/nesting_paper/analysis/outputs/tb_bigdel/plots/124321a0/minimap2/callsfilterpass_stats.tsv",sep="\t"))
# Num sequences mapped with edit distance of 0
#df_mapq30 %>% filter(NM == 0) %>% group_by(tool) %>% summarise(n()/765)
#df_lowmapq_minimap2 <- df_unfiltered_minimap2 %>% filter(!is.na(MAPQ) & MAPQ<=40)
#df_mapq_minimap2 <- df_unfiltered_minimap2
#df_mapq_minimap2$NM[df_mapq_minimap2$MAPQ < 30] <- NA
#max_NM <- max(df_mapq_minimap2$NM, na.rm = TRUE)
#plotted_dataset <-df_mapq_minimap2 %>% replace_na(list(NM = max_NM))
#ecdf_plot <- ggplot(plotted_dataset, aes(x = NM, colour = tool)) + stat_ecdf(geom="step", pad = FALSE)
#df_mapq_minimap2 %>% group_by(tool) %>% summarise(dist=mean(NM,na.rm=TRUE))
#long_missing_minimap2 <- df_mapq_minimap2 %>%
#  mutate(missing = as.numeric(is.na(NM))) %>%
#  select(sample, gene, tool, missing) %>%
#  spread(tool, missing)
#
#all_NAs_minimap2 <- df_mapq_minimap2 %>% filter(is.na(NM))
#gramtools_unmapped_minimap2 <- long_missing_minimap2 %>% filter(gramtools_df3b1583 == 1)
#graphtyper_unmapped_minimap2 <- long_missing_minimap2 %>% filter(graphtyper2 == 1)
#
#long_NM_minimap2 <- df_mapq_minimap2 %>%
#  select(sample, gene, tool, NM) %>%
#  spread(tool, NM)
#
#long_NM_minimap2 <- long_NM_minimap2 %>% mutate(gramtools_diff=gramtools_df3b1583 - graphtyper2)
#worse_than_gtyper2 <- long_NM_minimap2 %>% filter(gramtools_diff > 0)

