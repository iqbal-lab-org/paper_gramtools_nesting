library(dplyr)
library(purrr)
library(ggplot2)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot genotyping performance")
p <- add_argument(p, "input_tsv", help="")
p <- add_argument(p, "output_dir", help="")
argv <- parse_args(p)

get_classif <- function(has_call, is_correct){
 if (has_call){
    if (is_correct) return("TP")
    else return("FP")
  }
  else {
    if (is_correct) return("TN")
    else return("FN")
  }
}

get_recall <- function(classif){
  TP = sum(classif == "TP")
  TN = sum(classif == "TN")
  FN = sum(classif == "FN")
  return((TP + TN) / (TP + TN + FN))
}

get_precision <- function(classif){
  TP = sum(classif == "TP")
  TN = sum(classif == "TN")
  FP = sum(classif == "FP")
  return((TP + TN) / (TP + TN + FP))
}

data <- read_tsv(argv$input_tsv)
gmtools_commit <- basename(argv$output_dir)

# Add classification
data <- data %>% mutate(classif = pmap_chr(list(res_has_call, res_is_correct), get_classif))


# Compute and plot precision/recall
recalls <- data %>% group_by(prg, err_rate, fcov) %>% summarise(metric = "recall", score = get_recall(classif))
precisions <- data %>% group_by(prg, err_rate, fcov) %>% summarise(metric = "precision", score = get_precision(classif))
total <- rbind(recalls, precisions)
avg_recall = get_recall(data$classif)
avg_precision = get_precision(data$classif)
plot_title <- sprintf("avg recall=%.2f avg precision=%.2f, gmtools_commit: %s",avg_recall,avg_precision, gmtools_commit)

t <- ggplot(total, aes(prg,score)) + geom_bar(aes(fill=metric), stat="identity",position="dodge")
t <- t + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both) + labs(title = plot_title)
ggsave(file.path(argv$output_dir,"precision_recall.pdf"),width = 10, height = 8, plot=t)

# Plot correctness vs metric distributions
correctness <- data %>% filter(classif == "TP" | classif == "FP")
correctness <- correctness %>% mutate(log_GC = log(GC))

plot_title <- sprintf("genotype confidence x correctness, gmtools_commit: %s", gmtools_commit)
GC_boxplot <- ggplot(correctness, aes(classif, log_GC)) + geom_boxplot() + labs(title=plot_title)
GC_boxplot <- GC_boxplot + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
ggsave(file.path(argv$output_dir,"GC_distrib.pdf"),width = 8, height = 6, plot=GC_boxplot)


plot_title <- sprintf("genotype confidence percentile  x correctness, gmtools_commit: %s", gmtools_commit)
GCP_boxplot <- ggplot(correctness, aes(classif, GCP)) + geom_boxplot() + labs(title=plot_title)
GCP_boxplot <- GCP_boxplot + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
ggsave(file.path(argv$output_dir,"GCP_distrib.pdf"),width = 8, height = 6, plot=GCP_boxplot)
