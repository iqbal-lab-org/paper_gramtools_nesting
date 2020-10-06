library(ggplot2)
library(tibble)
library(readr)
library(dplyr)
library(argparser, quietly=TRUE)

p <- arg_parser("Plot prg_closest results")
p <- add_argument(p, "roc_tsv", help="")
p <- add_argument(p, "outdir", help="")
p <- add_argument(p, "condition", help="")

argv <- parse_args(p)

df = read_tsv(argv$roc_tsv)
#df = read_tsv("/home/brice/Desktop/main_PhD/analyses/nesting_paper/test/ROC_stats.tsv")
df = df %>% mutate(condition="GCP filter")
df2 = tibble(fpr=c(0, max(df$fpr)), tpr=c(0,max(df$fpr)), condition=c("y=x","y=x"))

df2 = bind_rows(df, df2)
roc_curve = ggplot(df2, aes(x=fpr,y=tpr, colour=condition)) + geom_line() + xlab("False positive rate") + 
  ylab("True positive rate") + theme(text=element_text(size = 13))
  fname <- sprintf("ROC_curve_%s.pdf", argv$condition)
ggsave(file.path(argv$outdir, fname), width=9, height=7, plot=roc_curve)

prec_recall_curve = ggplot(df, aes(y=one_minus_precision,x=tpr)) + geom_line() + ylab("1 - call precision") + 
  xlab("True positive rate") + theme(text=element_text(size = 13))
  fname <- sprintf("precision_recall_curve_%s.pdf", argv$condition)
ggsave(file.path(argv$outdir, fname), width=9, height=7, plot=prec_recall_curve)