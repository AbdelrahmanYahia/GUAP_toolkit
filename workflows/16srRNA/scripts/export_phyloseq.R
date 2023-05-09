library("optparse", quietly = T)
defaultW <- getOption("warn") 
options(warn = -1) 
############################################
option_list = list(
  make_option(c("-s", "--sv"), type="character", default=NULL, 
              help="input sequence table (table.qza)", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, 
              help="input metadata (sample-metadata.tsv)", metavar="character"),
  make_option(c("-x", "--taxonomy"), type="character", default=NULL, 
              help="input taxonomy (taxonomy.qza)", metavar="character"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="input tree (rooted-tree.qza)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="ps.RDS", 
              help="output file [default= %default]", metavar="character"),
  make_option(c("-w", "--working_dir"), type="character",default=getwd(), metavar="character")
);
#############################################
opt_parser = OptionParser(option_list=option_list,
);

opt = parse_args(opt_parser);

if (is.null(opt$sv)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [sequence table] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1)
}
if (is.null(opt$metadata)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [metadata] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1, runLast = TRUE)
}
if (is.null(opt$taxonomy)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [taxaonomy] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1, runLast = TRUE)
}
if (is.null(opt$tree)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [tree] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1, runLast = TRUE)
}

#######################################
library("phyloseq", quietly = T)
library("qiime2R", quietly = T)
#########################################
setwd(opt$working_dir)
metadata<-read_q2metadata(opt$metadata)
samdf <- metadata
rownames(metadata) <- metadata$SampleID
tree <- read_qza(opt$tree)$data
ps <- qza_to_phyloseq(
  features=opt$sv,
  tree=opt$tree, opt$taxonomy,
  metadata = opt$metadata
)
ps <- phyloseq(ps@otu_table, ps@tax_table, sample_data(metadata), read_tree(tree))
write(paste0("\033[0;", 32, "m","exporting RDS...","\033[0m"), stderr())
saveRDS(ps, opt$output)