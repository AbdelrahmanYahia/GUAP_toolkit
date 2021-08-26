############################################
#                                          #
#####       SCRIPT GUAP DE v1.1      #######
#####    developed by Abdelrhman    #######
#                                          #
############################################
# load libraries
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PCAtools")
library("WGCNA")
library("EnhancedVolcano")
library("tidyverse")
library("hrbrthemes")
library("viridis")
library("optparse")
###########  DE deseq2  ####################
DE_deseq2 <- function(target, pheno.data, cond1, cond2, design){
  # get normalized data for visualization
  dds = DESeqDataSetFromMatrix(countData = target, colData = pheno.data, design = ~Status)
  dds <- estimateSizeFactors(dds)
  target_normalized <- counts(dds, normalized=TRUE)
  target_normalized <- apply (target_normalized, 2,as.integer)
  rownames(target_normalized) <- rownames(target)
  as.data.frame(target_normalized) -> target_normalized
  # Run deseq2
  dds.run <- DESeq(dds)
  res <- results(dds.run, alpha=0.05)
  res = res[complete.cases(res),]
  summary(res)
  res.df <- as.data.frame(res)
  res.degs <- res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2),]
  outer <- list("dds" = dds ,"dds.run" = dds.run , "result.RAW" = res.df, "results.filterred" = res.degs, "normailzed.data" = target_normalized, "raw.deseq2.res" = res)
  return(outer)
  
}
########  Data preprocessing    ############
preproces <- function(all_samples_Count, clinical, c_zeros=1, s_zeros=1, S="s", C="c" ){
  # assign samples
  c.samples= clinical[ clinical$Status %in% c(C),]$Name
  s.samples= clinical[ clinical$Status %in% c(S),]$Name
  pheno.sub= clinical[clinical$Name %in% c(c.samples,s.samples), c("Name", "Status")]
  # get counts for samples 
  c.counts <- all_samples_Count[,c.samples]
  c.counts <- c.counts[apply(c.counts == 0, 1, sum) <= c_zeros, ]
  s.counts <- all_samples_Count[,s.samples]
  s.counts <- s.counts[apply(s.counts == 0, 1, sum) <= s_zeros, ]
  miRNA_names <- unique(c(rownames(c.counts), rownames(s.counts)))
  target <- all_samples_Count[,c(c.samples,s.samples)]
  target <- target[rownames(target) %in% miRNA_names, ]
  log_target <- log(target,10)
  log_target <- mutate_all(log_target, function(x) as.numeric(x))
  
  print(log_target %>% 
    pivot_longer(cols = names(log_target), 
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
    ggtitle("expression (log10)") +
    xlab("") +
    ylab(""))
  outer <- list("target" = target, "pheno.sub" = pheno.sub)
  return(outer)
}
##########  main run function  #############
RUN <- function(preproces){
  preproces <- preproces
  target <- preproces$target
  pheno.sub <- preproces$pheno.sub
  results_1 = DE_deseq2(target = target, pheno.data = pheno.sub, cond1 = C, cond2 = S, design = "Status")
  
  res.df <- as.data.frame(results_1$raw.deseq2.res)
  res.degs <- res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2),]
  exp.degs <- target[rownames(target) %in% rownames(res.degs), ]
  exp.degs <- exp.degs[ , order(names(exp.degs))]
  exp.degs <- exp.degs[order(row.names(exp.degs)) ,]
  outer <- list("exp.degs" = exp.degs, "res.df" = res.df, "deseq2_out" = results_1, "preproces" = preproces)
  return(outer)
}
########   command line run   ##############
cl_run <- function(all_human_Count, clinical, c_zeros=1, s_zeros=1, S="s", C="c", NM = "Sample"){
  NPOST <- preproces(all_human_Count, clinical,c_zeros,s_zeros, S=S, C=C)
  N <- RUN(NPOST)
  df1 <- data.frame( 
    meens = colMeans(NPOST$target),
    sums = colSums(NPOST$target)
  )
  df1 <- df1[order(rownames(df1)),]
  write.csv2(df1, file = paste(NM, "_Meanexpresion.csv", sep=""), quote = F )
  write.csv2(NPOST$target, file = paste(NM, "_raw-counts.csv", sep=""), quote = F )
  write.csv2(N$exp.degs, file = paste(NM, "_DEMs.csv", sep=""), quote = F )
  write.csv2(N$res.df, file = paste(NM, "_results.csv", sep=""), quote = F )
  get_RUN_figs(N)
}
###########  PCA function      #############
printpca <- function(preproces){
  Colkey <- data.frame(
    col <- c('red', 'forestgreen', "purple", 'black')
  )
  metad <- preproces$pheno.sub
  rownames(metad) <- metad$Name
  metad <- metad[ order(rownames(metad)),]
  target2 <- preproces$target
  target2 <- preproces$target[, order(names(preproces$target))]
  p <- pca(target2, metadata = metad)
  plotp <- biplot(p,
                  colby = 'Status', 
                  colkey = Colkey$col,
                  colLegendTitle = 'Sample',
                  legendPosition = 'top',
                  legendLabSize = 13, 
                  legendIconSize = 5.0)
  print(plotp)
}
#######   functions to get figures   #######
get_RUN_figs <- function(RUN){
  preproces <- RUN$preproces
  res <- RUN$deseq2_out$raw.deseq2.res
  res.degs <- RUN$deseq2_out$results.filterred
  exp.degs <- RUN$exp.degs
  
  printpca(preproces)
  
  plotClusterTreeSamples(datExpr=t(preproces$target))
  
  m2=scale(t(exp.degs),center=T,scale=T)
  m2=t(m2)
  pheatmap(m2, fontsize_row = 5,fontsize_col = 9,
           cluster_rows=T, cluster_cols=T, scale = "none",
           color = colorRampPalette(c("mintcream",  "steelblue4"))(50), border_color=NA)
  
  volc <- EnhancedVolcano(res,
                          title = "",
                          subtitle = bquote(italic(Volcanoplot)),
                          lab = rownames(res),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlab = bquote(~Log[2]~ 'fold change'),
                          pCutoff = 10e-4,
                          FCcutoff = 2.0,
                          pointSize = 3.0,
                          labSize = 3.5,
                          colAlpha = 1,
                          legendPosition = 'top',
                          legendLabSize = 6,
                          legendIconSize = 4.0,
                          drawConnectors = TRUE,
                          widthConnectors = 0.75)
  print(volc)
  
}
############################################
runfull<- function(all_human_Count, clinical,s=1, c=1, S="A", C="D", name= "AvD"){
  proc <- preproces(all_human_Count, clinical,s,c, S, C)
  RUNrez <- RUN(proc)
  get_RUN_figs(RUNrez, name)
  write.csv2(RUNrez$exp.degs, file = paste(name, ".csv", sep=""), quote = F )
}
############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input counts file", metavar="character"),
  make_option(c("-t", "--tab"), action="store_true",  default=FALSE,
              help="input is a tab separated file"),
  make_option(c("-c", "--clinical"), type="character", default=NULL, 
              help="input clincial file", metavar="character"),
  make_option(c("-S", "--sample"), type="character", default=NULL, 
              help="Sample name in clinical", metavar="character"),
  make_option(c("-C", "--control"), type="character", default=NULL, 
              help="Control name in clinical", metavar="character"),
  make_option(c("-z", "--sample-zeros"), type="integer", default=1, 
              help="Allowed number of zeros in sample group [default= %default]",
              metavar="number"),
  make_option(c("-x", "--control-zeros"), type="integer", default=1, 
              help="Allowed number of zeros in control group [default= %default]",
              metavar="number"),
  make_option(c("-o", "--out-dir"), type="character", default="out", 
              help="output file Directory [default= %default]", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="GUAP-DE-out", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-F", "--Full-run"), action="store_true", default=TRUE,
              help="Performe Full DE run [default]"),
  make_option(c("-Q", "--QC"), action="store_false", 
              dest="Full-run", help="Perfom QC only")
); 
#############################################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("No [input] was supplied", call.=FALSE)
}
if (is.null(opt$clinical)){
  print_help(opt_parser)
  stop("No [clinical] was supplied", call.=FALSE)
}
if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("No [sample] was supplied", call.=FALSE)
}
if (is.null(opt$control)){
  print_help(opt_parser)
  stop("No [control] was supplied", call.=FALSE)
}
############################################
# run analysis
re <- "_S([0-9]+)_L([0-9]+)_R(1|2)_([0-9]+).sam"
dirstr="OUT_miRNA_human_mapped_reads_"

clinical <- read.csv(opt$clinical)
names(clinical) <- c("Name", "Status") 
clinical$Name <- gsub(x = clinical$Name, pattern = "-", replacement = "_")
if ( opt$tab){
  all_human_Count <- read.csv(opt$input, row.names=1, sep = "\t")
}else{
  all_human_Count <- read.csv(opt$input, row.names=1, sep = ",")
}

names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = "\\.", replacement = "_")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = re, replacement = "")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = dirstr, replacement = "")


if ( opt$`Full-run` ) { 
  write("Performing Full DE Run...\n", stderr()) 

  proc <- preproces(all_human_Count, clinical,opt$`sample-zeros`,opt$`control-zeros`, S=opt$sample, C=opt$control)
  proc$pheno.sub
  RUNrez <- RUN(proc)
  get_RUN_figs(RUNrez)
  
} else {
  write("Performing QC only...\n", stderr())
  proc <- preproces(all_human_Count, clinical,opt$`sample-zeros`,opt$`control-zeros`, S=opt$sample, C=opt$control)
  printpca(proc)
  plotClusterTreeSamples(datExpr=t(proc$target))
}

