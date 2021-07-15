############################################
#                                          #
#####       SCRIPT GUAP DE v1.0      #######
#####    developped by Abdelrhman    #######
#                                          #
############################################
# load libraries
library("DESeq2")
library("edgeR")
library("DEGreport")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PCAtools")
library("WGCNA")

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
  plot(hclust((dist(t(preproces$target)))))
  plotClusterTreeSamples(datExpr=t(preproces$target))

  # Make a basic volcano plot
  par(mfrow=c(1,1))
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res.degs, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  
  m2=scale(t(exp.degs),center=T,scale=T)
  m2=t(m2)
  with(subset(res.degs, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  pheatmap(m2, fontsize_row = 5,fontsize_col = 9,
           cluster_rows=F, cluster_cols=F, scale = "none")
  
}
############################################
# get data
re <- "_S([0-9]+)_L([0-9]+)_R(1|2)_([0-9]+).sam"
dirstr="OUT_miRNA_human_mapped_reads_"
clinical <- read.csv("clinical.csv")
clinical$Name <- gsub(x = clinical$Name, pattern = "-", replacement = "_")

all_human_Count <- read.csv("final_human_counts.txt", row.names=1, sep = "\t")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = "\\.", replacement = "_")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = re, replacement = "")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = dirstr, replacement = "")

############################################
# direct_run_from_command_line <- cl_run(all_human_Count, clinical,1,3, S="ns", C="nc","Human_N")
############################################
# run analysis
proc <- preproces(all_human_Count, clinical,1,3, S="ns", C="nc")
RUNrez <- RUN(proc)
get_RUN_figs(RUNrez)

