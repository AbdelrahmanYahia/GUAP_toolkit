library("optparse",quietly = T)
library("dada2",quietly = T)
library("ggplot2",quietly = T)
library("viridis",quietly = T)
suppressMessages(library("hrbrthemes",quietly = T))
suppressMessages(library("tidyverse",quietly = T))
#############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input path", metavar="character"),

  make_option(c("-n", "--name"), type="character", default="/GUAP-16s-dada2", 
              help="Name of analysis [default= %default]", metavar="character"),

  make_option(c("-o", "--out-dir"), type="character", default="dada2", 
              help="output file Directory [default= %default]", metavar="character"),

  make_option(c("-t", "--trunclenf"), type="integer", default=0, 
              help="truncLen F ", metavar="number"),

  make_option(c("-T", "--trunclenr"), type="integer", default=0, 
              help="truncLen R ", metavar="number"),

  make_option(c("--metadata"), type="character", default=NULL, 
              help="input metadata (sample-metadata.tsv)", metavar="character"),

  make_option("--r1-pattern", type="character",
              help="DADA2 input files pattaren", metavar="character"),

  make_option("--r2-pattern", type="character",
              help="DADA2 input files pattaren", metavar="character"),

  make_option(c("-l", "--triml"), type="integer", default=0, 
              help="trim left ", metavar="number"),
  make_option(c("-r", "--trimr"), type="integer", default=0, 
              help="trim right ", metavar="number"),
  make_option(c("-e", "--maxeef"), type="integer", default=8, 
              help="maxEE F ", metavar="number"),
  make_option(c("-E", "--maxeer"), type="integer", default=9, 
              help="maxEE R ", metavar="number"),
  make_option(c("-m", "--minoverlap"), type="integer", default=10, 
              help="minoverlab to merge", metavar="number"),
  
  make_option(c("--chimethod"), type="character", default='consensus', 
              help="chimera method", metavar="character"),

  make_option(c("-p", "--threads"), type="integer", default=8, 
              help="Number of threads",
              metavar="number"),

  make_option(c("-a", "--assign-taxa"), action="store_true",
              help="Use dada2 to perform taxa identification", default=FALSE),
  
  make_option(c("--use-exsisting"), action="store_true",default=FALSE,
              help="Use dada2 to perform taxa identification"),
  
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
path <- opt$input 
outdir <- opt$`out-dir`
dir.create(paste0(outdir), showWarnings = FALSE)
setwd(outdir)

assigntaxa <- function(){
  if (opt$`assign-taxa`) {
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
  }
}

dir.create(paste0(outdir), showWarnings = FALSE)

fnFs <- sort(list.files(path, pattern=opt$`r1-pattern`, full.names = TRUE))
fnRs <- sort(list.files(path, pattern=opt$`r2-pattern`, full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(outdir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(outdir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","out",".RDS"))){
    out <- readRDS(paste0(outdir,"/","out",".RDS"))
  }else{
  write(paste0("\033[0;", 32, "m","Filter and Trim started...","\033[0m"), stderr())
  # write(paste0("\033[0;", 32, "m","parameters are: ",opt$trunclenf,opt$trunclenr,opt$maxeef,opt$maxeer,"\033[0m"), stderr())
  
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      maxN=0, truncQ=2, rm.phix=TRUE, 
                      maxEE=c(opt$maxeef,opt$maxeer),
                      truncLen=c(opt$trunclenf,opt$trunclenr),
                      trimLeft=opt$triml, trimRight=opt$trimr,
                      compress=TRUE, multithread=opt$threads, verbose=opt$verbose)
  saveRDS(out,paste0(outdir,"/","out",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","Filter and Trim started...","\033[0m"), stderr())
  # write(paste0("\033[0;", 32, "m","parameters are: ",opt$trunclenf,opt$trunclenr,opt$maxeef,opt$maxeer,"\033[0m"), stderr())
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                      maxN=0, truncQ=2, rm.phix=TRUE, 
                      maxEE=c(opt$maxeef,opt$maxeer),
                      truncLen=c(opt$trunclenf,opt$trunclenr),
                      trimLeft=opt$triml, trimRight=opt$trimr,
                      compress=TRUE, multithread=opt$threads, verbose=opt$verbose)
  saveRDS(out,paste0(outdir,"/","out",".RDS"))
}

if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","errF",".RDS"))){
    errF <- readRDS(paste0(outdir,"/","errF",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","Learn errors F started...","\033[0m"), stderr())
    errF <- suppressMessages(learnErrors(filtFs, multithread=opt$threads, nbases=90000, verbose=opt$verbose))
    saveRDS(errF,paste0(outdir,"/","errF",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","Learn errors F started...","\033[0m"), stderr())
  errF <- suppressMessages(learnErrors(filtFs, multithread=opt$threads, nbases=90000, verbose=opt$verbose))
  saveRDS(errF,paste0(outdir,"/","errF",".RDS"))
}

if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","errR",".RDS"))){
    errR <- readRDS(paste0(outdir,"/","errR",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","Learn errors R started...","\033[0m"), stderr())
    errR <- suppressMessages(learnErrors(filtRs, multithread=opt$threads, nbases=90000, verbose=opt$verbose))
    saveRDS(errR,paste0(outdir,"/","errR",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","Learn errors R started...","\033[0m"), stderr())
  errR <- suppressMessages(learnErrors(filtRs, multithread=opt$threads, nbases=90000, verbose=opt$verbose))
  saveRDS(errR,paste0(outdir,"/","errR",".RDS"))
}

if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","dadaFs",".RDS"))){
    dadaFs <- readRDS(paste0(outdir,"/","dadaFs",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","DADA infer F started...","\033[0m"), stderr())
    dadaFs <- dada(filtFs, err=errF, multithread=opt$threads, verbose=opt$verbose)
    saveRDS(dadaFs,paste0(outdir,"/","dadaFs",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","DADA infer F started...","\033[0m"), stderr())
  dadaFs <- dada(filtFs, err=errF, multithread=opt$threads, verbose=opt$verbose)
  saveRDS(dadaFs,paste0(outdir,"/","dadaFs",".RDS"))
}


if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","dadaRs",".RDS"))){
    dadaRs <- readRDS(paste0(outdir,"/","dadaRs",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","DADA infer R started...","\033[0m"), stderr())
    dadaRs <- dada(filtRs, err=errR, multithread=opt$threads, verbose=opt$verbose)
    saveRDS(dadaRs,paste0(outdir,"/","dadaRs",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","DADA infer R started...","\033[0m"), stderr())
  dadaRs <- dada(filtRs, err=errR, multithread=opt$threads, verbose=opt$verbose)
  saveRDS(dadaRs,paste0(outdir,"/","dadaRs",".RDS"))
}


if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","mergers",".RDS"))){
    mergers <- readRDS(paste0(outdir,"/","mergers",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","Merging started...","\033[0m"), stderr())
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,minOverlap=opt$minoverlap ,verbose=opt$verbose)
    saveRDS(mergers,paste0(outdir,"/","mergers",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","Merging started...","\033[0m"), stderr())
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,minOverlap=opt$minoverlap ,verbose=opt$verbose)
  saveRDS(mergers,paste0(outdir,"/","mergers",".RDS"))
}

    
if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","seqtab",".RDS"))){
    seqtab <- readRDS(paste0(outdir,"/","seqtab",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","ASV generation started...","\033[0m"), stderr())
    seqtab <- makeSequenceTable(mergers)
    saveRDS(seqtab,paste0(outdir,"/","seqtab",".RDS"))
  }
}else{
  write(paste0("\033[0;", 32, "m","ASV generation started...","\033[0m"), stderr())
  seqtab <- makeSequenceTable(mergers)
  saveRDS(seqtab,paste0(outdir,"/","seqtab",".RDS"))
}
    

if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","seqtab.nochim",".RDS"))){
    seqtab.nochim <- readRDS(paste0(outdir,"/","seqtab.nochim",".RDS"))
  }else{
    write(paste0("\033[0;", 32, "m","remove chimera started...","\033[0m"), stderr())
    seqtab.nochim <- removeBimeraDenovo(seqtab, method=opt$chimethod, multithread=opt$threads, verbose=opt$verbose)
    saveRDS(seqtab.nochim,paste0(outdir,"/","seqtab.nochim",".RDS"))
    
    chim_perc <- (sum(seqtab.nochim)/sum(as.data.frame(out)$reads.in)) * 100
    write(paste0("\033[0;", 33, "m","filterred reads % = ",round(chim_perc,2),"\033[0m"), stderr())
  }
}else{
  write(paste0("\033[0;", 32, "m","remove chimera started...","\033[0m"), stderr())
  seqtab.nochim <- removeBimeraDenovo(seqtab, method=opt$chimethod, multithread=opt$threads, verbose=opt$verbose)
  saveRDS(seqtab.nochim,paste0(outdir,"/","seqtab.nochim",".RDS"))
  
  chim_perc <- (sum(seqtab.nochim)/sum(as.data.frame(out)$reads.in)) * 100
  print(paste("filterred reads % = ",round(chim_perc,2)))
  write(paste0("\033[0;", 33, "m","filterred reads % = ",round(chim_perc,2),"\033[0m"), stderr())
}


getN <- function(x) sum(getUniques(x))

if (length(sample.names) > 1 ){
  track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
} else {
  track <- cbind(out, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
}

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
write(paste0("\033[0;", 32, "m","writing output...","\033[0m"), stderr())

write.csv(track, file = paste0(outdir,"/stats.csv"), quote = FALSE)
write.table(t(seqtab.nochim), paste0(outdir,"/seqtab-nochim.txt"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, fout=paste0(outdir,'/rep-seqs.fna'), ids=colnames(seqtab.nochim))


if (opt$`use-exsisting`){
  if ( file.exists(paste0(outdir,"/","taxa",".RDS"))){
    taxa <- readRDS(paste0(outdir,"/","taxa",".RDS"))
  }else{
    assigntaxa()
  }
}else{
  assigntaxa()
}

################       downstream analysis    ##################################
write(paste0("\033[0;", 32, "m","Generating figures...","\033[0m"), stderr())
t.track <- t(track)
t.track <- as.data.frame(t.track)
tryCatch({
    tbox <- (t.track) %>% 
      pivot_longer(cols = names((t.track)), 
                    names_to = "gene", 
                    values_to = "expression") %>% 
      ggplot( aes(x=gene, y=expression, fill=gene)) +
      geom_boxplot() +
      scale_fill_viridis(discrete = TRUE, alpha=0.6) +
      theme_ipsum() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      ggtitle("Counts") +
      xlab("") +
      ylab("")

    ggsave(paste0(outdir,"/",opt$name, "-stats_box.png"), tbox ,device = png, width = 30, height = 10, units = "in" ,dpi = 350)
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error in Plotting Box plot for DADA2 stats: ","\033[0m",err), stderr())
  }
)

tryCatch({
    t.track2 <- as.data.frame(t(data.frame(
      input=track$input,
      filttered=track$nonchim
    )))

    names(t.track2) <- names(t.track)
    p11 <- (t.track2) %>% 
      pivot_longer(cols = names((t.track)), 
                    names_to = "gene", 
                    values_to = "expression")

    ns <- length(sample.names)

    p11$data <- c(rep(c("input"), times = ns),rep(c("output"), times = ns))

    t2box <- ggplot(p11, aes(fill=data, y=expression, x=gene)) + 
      geom_bar(position="stack", stat="identity") +
      scale_fill_viridis(discrete = T) +
      ggtitle("Counts filttered") +
      theme_ipsum() +
      xlab("")

    ggsave(paste0(outdir,"/",opt$name, "-comp_box.png"), t2box ,device = png, width = 25, height = 10, units = "in" ,dpi = 350)
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error in Plotting bar plot for DADA2 stats: ","\033[0m",err), stderr())
  }
)

write(paste0("\033[0;", 32, "m","Done...","\033[0m"), stderr())