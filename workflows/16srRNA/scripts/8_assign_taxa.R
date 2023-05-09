library("optparse",quietly = T)
library("dada2",quietly = T)

#############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input path", metavar="character"),
  make_option(c("-o", "--out-dir"), type="character", default="dada2", 
              help="output file Directory [default= %default]", metavar="character"),
  
  make_option(c("-n", "--name"), type="character", default="/GUAP-16s-dada2", 
              help="Name of analysis [default= %default]", metavar="character"),

  make_option(c("-p", "--threads"), type="integer", default=8, 
              help="Number of threads",
              metavar="number"),
  make_option(c("-c", "--classifier"),type="character",  metavar="character",
              help="Use dada2 to perform taxa identification"),
              
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Verbose [default]"),
  make_option(c("-s", "--silent"), action="store_false", 
              dest="verbose", help="verbose = FALSE")
); 
#############################################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No [input] was supplied", call.=FALSE)
}

############################################
seqtab.nochim <- readRDS(paste0(outdir,"/","seqtab.nochim",".RDS"))

write(paste0("\033[0;", 32, "m","Assigning taxonomy...","\033[0m"), stderr())
taxa <- assignTaxonomy(seqtab.nochim, opt$classifier, multithread=opt$threads)
saveRDS(taxa,paste0(outdir,"/","taxa",".RDS"))
write(paste0("\033[0;", 32, "m","Exporting taxa...","\033[0m"), stderr())
library("phyloseq")
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                tax_table(taxa))
tax<-as(tax_table(ps),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL
write.table(tax, paste0(outdir,"/","tax.txt"), quote=FALSE, col.names=FALSE, sep="\t")
