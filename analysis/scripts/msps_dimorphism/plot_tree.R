library(argparser, quietly=TRUE)

p <- arg_parser("Plot phylo tree")
p <- add_argument(p, "phylo_tree", help="in Newick format")
p <- add_argument(p, "output_fname", help="will save with extension of .pdf")

argv <- parse_args(p)

require(ape)
require(ggtree)

input_tree <- read.tree(argv$phylo_tree)
display <- ggtree(input_tree) + ggtitle(paste("Tree: ", basename(argv$phylo_tree)))
ggsave(paste(argv$output_fname, ".pdf", sep=""), display)