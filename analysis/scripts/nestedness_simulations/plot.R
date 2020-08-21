#Plots precision/recall of calls and genotype confidence/genotype confidence percentile distributions

library(ggplot2)
library(tibble)
library(dplyr)
library(argparser, quietly=TRUE)
library(rjson)

p <- arg_parser("Plot genotyping performance")
p <- add_argument(p, "input_tsv", help="")
p <- add_argument(p, "output_dir", help="")
p <- add_argument(p, "dataset_name", help="")
p <- add_argument(p, "gramtools_commit", help="")

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
plot_prec_recall <- function(data, title, argv, prg_name, gmtools_commit, faceted = TRUE){
  recalls <- data %>% group_by(nesting, err_rate, fcov) %>% summarise(metric = "recall", score = get_recall(classif))
  precisions <- data %>% group_by(nesting, err_rate, fcov) %>% summarise(metric = "precision", score = get_precision(classif))
  total <- rbind(recalls, precisions)
  avg_recall = get_recall(data$classif)
  avg_precision = get_precision(data$classif)
  #plot_title <- sprintf("dataset: %s, gramtools_commit: %s", prg_name, gmtools_commit)
  plot_title <- ""
  plot_subtitle <-sprintf("average recall=%.3f average precision=%.3f", avg_recall,avg_precision)
  
  prec_recall <- ggplot(total, aes(metric,score)) + geom_bar(aes(fill=nesting), stat="identity",position="dodge") + 
    labs(title = plot_title, subtitle=plot_subtitle)
  if (faceted){
    prec_recall <- prec_recall + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
    title<-paste(argv$dataset_name,title,sep="_")
  }
  ggsave(file.path(argv$output_dir,title),width = 8, height = 6, plot=prec_recall)
}

# Plot correctness vs metric distributions
plot_GC <- function(data, title, argv, prg_name, gmtools_commit, faceted = TRUE, log_scale = FALSE){
  correctness <- data %>% filter(classif == "TP" | classif == "FP")
  if (log_scale)
    correctness <- correctness %>% mutate(GC = log(GC))
  
  #plot_title <- sprintf("dataset: %s, gramtools_commit: %s", prg_name, gmtools_commit)
  plot_title <- ""
  GC_boxplot <- ggplot(correctness, aes(classif, GC)) + geom_boxplot(aes(fill=nesting)) + 
    labs(title=plot_title) + xlab("Call classification") + ylab("Genotype confidence")
  if (faceted){
    GC_boxplot <- GC_boxplot + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
    title<-paste(argv$dataset_name,title,sep="_")
  }
    
  ggsave(file.path(argv$output_dir,title),width = 8, height = 6, plot=GC_boxplot)
}


plot_GCP <- function(data, title, argv, prg_name, gmtools_commit, faceted = TRUE){
  correctness <- data %>% filter(classif == "TP" | classif == "FP")
  correctness <- correctness %>% mutate(log_GC = log(GC))
  
  #plot_title <- sprintf("dataset: %s\ngenotype confidence percentile  x correctness, gramtools_commit: %s", prg_name, gmtools_commit)
  plot_title <- ""
  GCP_boxplot <- ggplot(correctness, aes(classif, GCP)) + geom_boxplot(aes(fill=nesting)) + labs(title=plot_title)
  if (faceted)
    GCP_boxplot <- GCP_boxplot + facet_grid(cols=vars(err_rate), rows=vars(fcov), labeller = label_both)
  
  title<-paste(argv$dataset_name,title,sep="_")
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
gmtools_commit <- argv$gramtools_commit
prg_name = argv$dataset_name


# Add classification
data <- data %>% mutate(classif = pmap_chr(list(res_has_call, res_is_correct), get_classif))

prg_data <- data %>% filter(prg == prg_name) 
print(sprintf("Num filtered out ambiguous sites for prg %s: %d out of %d",
              prg_name, table(prg_data$ambiguous)[2], dim(prg_data)[1]))

data_noambi <- filter(data, ambiguous == 0)
prg_data <- data_noambi %>% filter(prg == prg_name) 

data_lvl1 <- filter(prg_data, lvl_1 == 1)
all_data_fixed <- filter(data_noambi, err_rate==0)

plot_prec_recall(data_noambi, "precision_recall.pdf", argv, prg_name, gmtools_commit)
plot_prec_recall(data_lvl1, "precision_recall_lvl1.pdf", argv, prg_name, gmtools_commit)
plot_prec_recall(all_data_fixed, "precision_recall_all.pdf", argv,  "all genes", gmtools_commit, faceted = FALSE)

plot_GC(data_noambi, "GC_distrib.pdf", argv, prg_name, gmtools_commit)
plot_GC(all_data_fixed, "GC_distrib_all.pdf", argv, "all genes", gmtools_commit, faceted=FALSE)

plot_GCP(data_noambi, "GCP_distrib.pdf", argv, prg_name, gmtools_commit)
