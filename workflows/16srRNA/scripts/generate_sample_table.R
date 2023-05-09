library("optparse",quietly = T)
library("tools",quietly = T)
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input path", metavar="character"),
  make_option("--r1-pattern", type="character",
              help="DADA2 input files pattaren", metavar="character"),
  make_option("--r2-pattern", type="character",
              help="DADA2 input files pattaren", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output path", metavar="character")
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
path <- opt$input 
path_out <- opt$output

fnFs <- sort(list.files(file_path_as_absolute(path), pattern=opt$`r1-pattern`, full.names = TRUE))
fnRs <- sort(list.files(file_path_as_absolute(path), pattern=opt$`r2-pattern`, full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
header <- c("Sample-id", "forward-absolute-filepath", "reverse-absolute-filepath")


df <- data.frame(
  `Sample-id` = sample.names,
  `forward-absolute-filepath` = fnFs,
  `reverse-absolute-filepath` = fnRs
)

names(df) <- header
write.table(df, paste0(path_out, "/", "samples.tsv"), sep = "\t", quote = F, row.names = F)
