############################################
library("optparse", quietly = T)
suppressMessages(library("tools", quietly = T))

defaultW <- getOption("warn") 
options(warn = -1) 

############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input counts file", metavar="character"),
  make_option(c("-t", "--tab"), action="store_true",  default=FALSE,
              help="input is a tab separated file"),
  make_option(c("-d", "--dirstr"), type="character", default="star/", 
              help="Dir str to remove [default= %default]", metavar="character"),
  make_option(c("-r", "--re"), type="character", default="/Aligned.sortedByCoord.out.bam", 
              help="RE pattern to remove [default= %default]", metavar="character")
); 

#############################################

opt_parser = OptionParser(option_list=option_list);

opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [INPUT] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1)
}

if ( opt$tab){
  all_human_Count <- read.csv(file_path_as_absolute(opt$input), row.names=1, sep = "\t", check.names = FALSE)
}else{
  all_human_Count <- read.csv(file_path_as_absolute(opt$input), row.names=1, sep = ",", check.names = FALSE)
}
re <- opt$re
dirstr <- opt$dirstr
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = "\\.", replacement = "_")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = re, replacement = "")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = dirstr, replacement = "")

write.table(all_human_Count, paste0(opt$input), quote=FALSE, col.names=TRUE, sep="\t", row.names = TRUE)
