library(tidyverse)
library(argparser, quietly=TRUE)
library(rjson)

p <- arg_parser("Plot genotyping performance")
p <- add_argument(p, "input_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "nested_json", help="")

#______Classifying calls____#
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


#_______Plots_______#
#Compute and plot precision/recall
plot_prec_recall <- function(data, title, argv, prg_name, gmtools_commit){
  recalls <- data %>% group_by(nesting, err_rate, fcov) %>% summarise(metric = "recall", score = get_recall(classif))
  precisions <- data %>% group_by(nesting, err_rate, fcov) %>% summarise(metric = "precision", score = get_precision(classif))
  total <- rbind(recalls, precisions)
  avg_recall = get_recall(data$classif)
  avg_precision = get_precision(data$classif)
  plot_title <- sprintf("dataset: %s \navg recall=%.2f avg precision=%.2f, gmtools_commit: %s",prg_name, avg_recall,avg_precision, gmtools_commit)
  
  prec_recall <- ggplot(total, aes(metric,score)) + geom_bar(aes(fill=nesting), stat="identity",position="dodge")
  prec_recall <- prec_recall + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both) + labs(title = plot_title)
  
  ggsave(file.path(argv$output_dir,title),width = 10, height = 8, plot=prec_recall)
}

# Plot correctness vs metric distributions
plot_GC <- function(data, title, argv, prg_name, gmtools_commit){
  correctness <- data %>% filter(classif == "TP" | classif == "FP")
  correctness <- correctness %>% mutate(log_GC = log(GC))
  
  plot_title <- sprintf("dataset: %s \ngenotype confidence x correctness, gmtools_commit: %s", prg_name, gmtools_commit)
  GC_boxplot <- ggplot(correctness, aes(classif, log_GC)) + geom_boxplot(aes(fill=nesting)) + labs(title=plot_title)
  GC_boxplot <- GC_boxplot + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
  ggsave(file.path(argv$output_dir,title),width = 8, height = 6, plot=GC_boxplot)
}


plot_GCP <- function(data, title, argv, prg_name, gmtools_commit){
  correctness <- data %>% filter(classif == "TP" | classif == "FP")
  correctness <- correctness %>% mutate(log_GC = log(GC))
  
  plot_title <- sprintf("dataset: %s\ngenotype confidence percentile  x correctness, gmtools_commit: %s", prg_name, gmtools_commit)
  GCP_boxplot <- ggplot(correctness, aes(classif, GCP)) + geom_boxplot(aes(fill=nesting)) + labs(title=plot_title)
  GCP_boxplot <- GCP_boxplot + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
  ggsave(file.path(argv$output_dir,title),width = 8, height = 6, plot=GCP_boxplot)
}

print_prec_recall <- function(data){
  data <- data %>% mutate(classif = pmap_chr(list(res_has_call, res_is_correct), get_classif))
  recalls <- data %>% group_by(nesting, err_rate, fcov) %>% summarise(metric = "recall", score = get_recall(classif))
  precisions <- data %>% group_by(nesting, err_rate, fcov) %>% summarise(metric = "precision", score = get_precision(classif))
  print(precisions)
  print(recalls)
}


#_____Main code: data load, process, plot____#
argv <- parse_args(p)
data <- read_tsv(argv$input_tsv)
gmtools_commit <- basename(argv$output_dir)
nested_json <- fromJSON(file=argv$nested_json)

prg_name = "DBLMSPs"
data <- data %>% filter(prg == prg_name) 

# Add classification
data <- data %>% mutate(classif = pmap_chr(list(res_has_call, res_is_correct), get_classif))

print(sprintf("Num filtered out ambiguous sites: %d out of %d",table(data$ambiguous)[2], dim(data)[1]))
data_noambi <- filter(data, ambiguous == 0)
data_lvl1 <- filter(data_noambi, lvl_1 == 1)

plot_prec_recall(data_noambi, "precision_recall.pdf", argv, prg_name, gmtools_commit)
plot_prec_recall(data_lvl1, "precision_recall_lvl1.pdf", argv, prg_name, gmtools_commit)

plot_GC(data_noambi, "GC_distrib.pdf", argv, prg_name, gmtools_commit)
plot_GC(data_lvl1, "GC_distrib_lvl1.pdf", argv, prg_name, gmtools_commit)

plot_GCP(data_noambi, "GCP_distrib.pdf", argv, prg_name, gmtools_commit)
plot_GCP(data_lvl1, "GCP_distrib_lvl1.pdf", argv, prg_name, gmtools_commit)
