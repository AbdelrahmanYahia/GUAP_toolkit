############################################
#                                          #
#####       SCRIPT GUAP DE v1.1      #######
#####    developed by Abdelrhman    #######
#                                          #
############################################
library("optparse", quietly = T)
defaultW <- getOption("warn") 
options(warn = -1) 
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

  make_option(c("-d", "--dirstr"), type="character", default="OUT_miRNA_human_mapped_reads_", 
              help="Dir str to remove [default= %default]", metavar="character"),
  make_option(c("-r", "--re"), type="character", default="_S([0-9]+)_L([0-9]+)_R(1|2)_([0-9]+).sam", 
              help="RE pattern to remove [default= %default]", metavar="character"),
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

opt_parser = OptionParser(option_list=option_list,
                          usage = "usage: GUAP DE -i [INPUT] -c [clinical] -S [sample] -C [control] -o [out-dir]",
                          prog = "DE",
                          description = "Performs DE on counts table",
);

opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [INPUT] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1)
}
if (is.null(opt$clinical)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [CLINICAL] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1, runLast = TRUE)
}
if (is.null(opt$sample)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [sample] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1, runLast = TRUE)
}
if (is.null(opt$control)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [control] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1, runLast = TRUE)
}

############################################
write(paste0("\033[0;", 33, "m","loading libraries...","\033[0m"), stderr())
suppressMessages(library("DESeq2",quietly = T))
suppressMessages(library("ggplot2", quietly = T))
suppressMessages(library("pheatmap", quietly = T))
suppressMessages(library("RColorBrewer", quietly = T))
suppressMessages(library("PCAtools", quietly = T))
suppressMessages(library("WGCNA", quietly = T))
suppressMessages(library("EnhancedVolcano", quietly = T))
suppressMessages(library("tidyverse", quietly = T))
suppressMessages(library("hrbrthemes", quietly = T))
suppressMessages(library("viridis", quietly = T))
suppressMessages(library("DEGreport", quietly = T))
suppressMessages(library("tools", quietly = T))
###########  DE deseq2  ####################
DE_deseq2 <- function(target, pheno.data, cond1, cond2){
  # get normalized data for visualization
  dds = suppressMessages(DESeqDataSetFromMatrix(countData = target, colData = pheno.data, design = ~Status))

  # Run deseq2
  dds.run <- suppressMessages(DESeq(dds))
  res <- results(dds.run, alpha=0.05,contrast = c("Status",cond2,cond1))
  res = res[complete.cases(res),]
  summary(res)
  res.df <- as.data.frame(res)
  res.degs <- res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2),]
  target_normalized <- counts(dds.run, normalized=TRUE)
  target_normalized <- apply (target_normalized, 2,as.integer)
  rownames(target_normalized) <- rownames(target)
  as.data.frame(target_normalized) -> target_normalized
  outer <- list("dds" = dds ,"dds.run" = dds.run , "result.RAW" = res.df, "results.filterred" = res.degs, "normailzed.data" = target_normalized, "raw.deseq2.res" = res)
  return(outer)
}
########  Data preprocessing    ############
preproces <- function(all_samples_Count, clinical, c_zeros=1, s_zeros=1, S="s", C="c",name ){
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
  
  Tplot <- suppressMessages(log_target %>% 
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
  suppressMessages(ggsave(paste0(name,"_Sample_counts.png"),Tplot))
  outer <- list("target" = target, "pheno.sub" = pheno.sub, "S"=S, "C"=C)
  return(outer)
}
##########  main run function  #############
RUN <- function(preproces,name){
  preproces <- preproces
  target <- preproces$target
  pheno.sub <- preproces$pheno.sub
  S=preproces$S
  C=preproces$C

  results_1 = DE_deseq2(target = target, pheno.data = pheno.sub, cond1 = C, cond2 = S)
  
  res.df <- as.data.frame(results_1$raw.deseq2.res)
  res.degs <- res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2),]
  exp.degs <- target[rownames(target) %in% rownames(res.degs), ]
  exp.degs <- exp.degs[ , order(names(exp.degs))]
  exp.degs <- exp.degs[order(row.names(exp.degs)) ,]
  if (nrow(exp.degs) > 1 ){
    write.csv(res.degs, paste0(name,"_DEMs.csv"), quote=F)
  }else{
    write(paste0("\033[0;", 31, "m","No DEMs Found...","\033[0m"), stderr())
  }
  write.csv(res.df, paste0(name,"_Results.csv"), quote=F)
  outer <- list("exp.degs" = exp.degs, "res.df" = res.df, "deseq2_out" = results_1, "preproces" = preproces)
  return(outer)
}
###########  PCA function      #############
printpca <- function(preproces, name){
  Colkey <- data.frame(
    col <- c('red', 'forestgreen', "purple", 'black')
  )
  metad <- preproces$pheno.sub
  rownames(metad) <- metad$Name
  metad <- metad[ order(rownames(metad)),]
  target2 <- preproces$target
  target2 <- preproces$target[, order(names(preproces$target))]
  p <- pca(target2, metadata = metad)
  write.csv(loadings(p), paste0(name,"pca_data.csv"), quote=F)

  plotp <- suppressMessages(biplot(p,
                  colby = 'Status', 
                  colkey = Colkey$col,
                  colLegendTitle = 'Sample',
                  legendPosition = 'top',
                  legendLabSize = 13, 
                  legendIconSize = 5.0))
  suppressMessages(print(plotp))
}
#######   functions to get figures   #######
get_RUN_figs <- function(RUN, name){
  preproces <- RUN$preproces
  res <- RUN$deseq2_out$raw.deseq2.res
  res.degs <- RUN$deseq2_out$results.filterred
  exp.degs <- RUN$exp.degs
  dds <- RUN$deseq2_out$dds.run
  svg(paste(name, "_PCA.svg", sep=""),  width=14, height=9)
  printpca(preproces, name)
  dev.off()
  svg(paste(name, "_Cluster_tree.svg", sep=""),  width=14, height=8)

  plotClusterTreeSamples(datExpr=t(preproces$target))
  dev.off()
  m2=scale(t(exp.degs),center=T,scale=T)
  m2=t(m2)
  write.csv(m2, paste0(name,"_HM_data.csv"), quote=F)

  if (nrow(exp.degs) > 1 ){
      svg(paste(name, "_Heatmap.svg", sep=""))
      hm <- pheatmap(m2, fontsize_row = 5,fontsize_col = 9,
           cluster_rows=F, cluster_cols=F, scale = "none",
           color = colorRampPalette(c("mintcream",  "steelblue4"))(50), border_color=NA)
      print(hm)
      dev.off()
      n = as.integer(nrow(exp.degs))
      svg(paste(name, "_demplot.svg", sep=""),  width=10, height=5)
      print(degPlot(dds = dds, res = res, n = n, xs = "Status", group="Status" ))
      dev.off()
  }

  svg(paste(name, "_volcanoplot_E_non_Filtered.svg", sep=""),  width=15, height=10)
  volc <- EnhancedVolcano(res,
                          title = "",
                          subtitle = bquote(italic(Volcanoplot)),
                          lab = rownames(res),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlab = bquote(~Log[2]~ 'fold change'),
                          pCutoff = 10e-4,
                          FCcutoff = 1.0,
                          pointSize = 3.0,
                          labSize = 3.5,
                          colAlpha = 1,
                          legendPosition = 'top',
                          legendLabSize = 6,
                          legendIconSize = 4.0,
                          drawConnectors = TRUE,
                          widthConnectors = 0.75)
  
  suppressMessages(print(volc))
  dev.off()

}
############################################
write(paste0("\033[0;", 33, "m","reading files...","\033[0m"), stderr())
# run analysis
re = opt$re 
dirstr = opt$dirstr
wd = opt$wd
current = getwd()
outdir <- opt$`out-dir`
dir.create(paste0(outdir), showWarnings = FALSE)

nameexp <- opt$name

clinical <- read.csv(file_path_as_absolute(opt$clinical), check.names = FALSE, sep = '\t')
names(clinical) <- c("Name", "Status") 


if ( opt$tab){
  all_human_Count <- read.csv(file_path_as_absolute(opt$input), row.names=1, sep = "\t", check.names = FALSE)
}else{
  all_human_Count <- read.csv(file_path_as_absolute(opt$input), row.names=1, sep = ",", check.names = FALSE)
}

names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = "\\.", replacement = "_")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = re, replacement = "")
names(all_human_Count) <- gsub(x = names(all_human_Count), pattern = dirstr, replacement = "")
setwd(outdir)

if ( opt$`Full-run` ) { 
  write(paste0("\033[0;", 33, "m","Preprocessing data...","\033[0m"), stderr())
  proc <- suppressMessages(preproces(all_human_Count, clinical,opt$`sample-zeros`,opt$`control-zeros`, S=opt$sample, C=opt$control,nameexp))
  proc$pheno.sub
  write(paste0("\033[0;", 33, "m","Performing DE...","\033[0m"), stderr())
  RUNrez <- suppressMessages(RUN(proc, nameexp))
  write(paste0("\033[0;", 33, "m","Generating figures...","\033[0m"), stderr())
  suppressMessages(get_RUN_figs(RUNrez, nameexp))
  write(paste0("\033[0;", 32, "m","DONE FULL RUN...","\033[0m"), stderr())

} else {
  write("Performing QC only...\n", stderr())
  proc <- suppressMessages(preproces(all_human_Count, clinical,opt$`sample-zeros`,opt$`control-zeros`, S=opt$sample, C=opt$control,nameexp))
  suppressMessages(printpca(proc, nameexp))
  suppressMessages(plotClusterTreeSamples(datExpr=t(proc$target)))
  RUNrez <- suppressMessages(RUN(proc, nameexp))
  proc <- suppressMessages(preproces(RUNrez$deseq2_out$normailzed.data, clinical,opt$`sample-zeros`,opt$`control-zeros`, S=opt$sample, C=opt$control,nameexp))
  suppressMessages(plotClusterTreeSamples(datExpr=t(proc$target)))
  suppressMessages(printpca(proc, nameexp))
  write(paste0("\033[0;", 32, "m","DONE QC...","\033[0m"), stderr())
}
options(warn = defaultW)
