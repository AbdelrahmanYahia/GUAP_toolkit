library("optparse",quietly = T)
library("dada2",quietly = T)
library("ggplot2",quietly = T)
library("viridis",quietly = T)
library("hrbrthemes",quietly = T)
suppressMessages(library("tidyverse",quietly = T))
#############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input path", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="/GUAP-16s-dada2", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-o", "--out-dir"), type="character", default="dada2", 
              help="output file Directory [default= %default]", metavar="character"),
  make_option(c("-t", "--trunclenf"), type="integer", default=290, 
              help="truncLen F ", metavar="number"),
  make_option(c("-T", "--trunclenr"), type="integer", default=220, 
              help="truncLen R ", metavar="number"),

  make_option(c("-e", "--maxeef"), type="integer", default=8, 
              help="maxEE F ", metavar="number"),
  make_option(c("-E", "--maxeer"), type="integer", default=9, 
              help="maxEE R ", metavar="number"),

  make_option(c("-p", "--threads"), type="integer", default=8, 
              help="Number of threads",
              metavar="number"),

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
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
write(paste0("\033[0;", 32, "m","Filter and Trim started...","\033[0m"), stderr())

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(opt$trunclenf,opt$trunclenr),
                     maxN=0, truncQ=2, rm.phix=TRUE, maxEE=c(opt$maxeef,opt$maxeer),
                     compress=TRUE, multithread=opt$threads, verbose=opt$verbose)

write(paste0("\033[0;", 32, "m","Learn errors started...","\033[0m"), stderr())
errF <- learnErrors(filtFs, multithread=opt$threads, nbases=90000, verbose=opt$verbose)
errR <- learnErrors(filtRs, multithread=opt$threads, nbases=90000, verbose=opt$verbose)

suppressMessages(plotErrors(errF, nominalQ=TRUE))

dadaFs <- dada(filtFs, err=errF, multithread=opt$threads, verbose=opt$verbose)
dadaRs <- dada(filtRs, err=errR, multithread=opt$threads, verbose=opt$verbose)

write(paste0("\033[0;", 32, "m","Merging started...","\033[0m"), stderr())
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=opt$verbose)
seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths

write(paste0("\033[0;", 32, "m","remove chimera started...","\033[0m"), stderr())
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=opt$threads, verbose=opt$verbose)

chim_perc <- (sum(seqtab.nochim)/sum(seqtab)) * 100
print(paste("chimera % = ",round(chim_perc,2)))
write(paste0("\033[0;", 33, "m","chimera % = ",round(chim_perc,2),"\033[0m"), stderr())

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- as.data.frame(track)
write(paste0("\033[0;", 32, "m","writing output...","\033[0m"), stderr())

write.csv(track, file = paste0(outdir,"/stats.csv"), quote = FALSE)
write.table(t(seqtab.nochim), paste0(outdir,"/seqtab-nochim.txt"), sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab.nochim, fout=paste0(outdir,'/rep-seqs.fna'), ids=colnames(seqtab.nochim))

save.image(paste0(outdir,"/",opt$name,".RData"))
################       downstream analysis    ##################################
###### check counts #####
write(paste0("\033[0;", 32, "m","Generating figures...","\033[0m"), stderr())
t.track <- t(track)
t.track <- as.data.frame(t.track)
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

################################################################################

write(paste0("\033[0;", 32, "m","Done...","\033[0m"), stderr())
