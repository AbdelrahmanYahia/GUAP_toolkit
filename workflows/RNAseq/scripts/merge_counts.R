library("optparse", quietly = T)
library(tibble, quietly = T)
defaultW <- getOption("warn") 
options(warn = -1) 
#############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output file Directory [default= %default]", metavar="character"),
  make_option(c("-s", "--srtandedness"),type="integer",  metavar="integer",default=1, 
              help="1: not stranded, 2: stranded, 3: reverse")
); 
#############################################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No [input] was supplied", call.=FALSE)
}
if (is.null(opt$output)){
  print_help(opt_parser)
  stop("No [output] was supplied", call.=FALSE)
}
############################################

input <- opt$input
dirs <- dir(input)
strand <- opt$srtandedness

merged_data <- data.frame()

for (sample in dirs) {
  filename <- sample
  path <- paste0(input, "/", filename, "/ReadsPerGene.out.tab")
  file1 <- read.csv(path, sep = "\t", header = FALSE, skip = 4, row.names = 1)
  names(file1) <- c("NON-strand", "Strand", "reverse")
  file_temp <- as.data.frame(file1[, strand])
  rownames(file_temp) <- rownames(file1)
  names(file_temp) <- filename
  
  if (nrow(merged_data) == 0) {
    merged_data <- file_temp
  } else {
    merged_data <- cbind(merged_data, file_temp)
  }
}

# move rownames to first column
merged_data <- rownames_to_column(as.data.frame(merged_data), var = "Genes")
write.table(merged_data, paste0(opt$output), quote=FALSE, col.names=TRUE, sep="\t", row.names = FALSE)