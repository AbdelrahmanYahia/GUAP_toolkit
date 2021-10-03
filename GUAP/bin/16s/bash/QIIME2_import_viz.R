library("optparse", quietly = T)
defaultW <- getOption("warn") 
options(warn = -1) 
#######################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input counts file", metavar="character"),
  make_option(c("-r", "--rare"), type="integer", default=2500, 
              help="input clincial file", metavar="number"),
  make_option(c("-n", "--name"), type="character", default="GUAP-DE-out", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-c", "--condition"), type="character", default="condition", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-S", "--sample"), type="character", default="sample", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-C", "--control"), type="character", default="control", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-w", "--working_dir"), type="character",default=getwd(), metavar="character")
);
#######################################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
write(paste0("\033[0;", 32, "m","Generating Figures started...","\033[0m"), stderr())
#######################################
suppressMessages(library("phyloseq", quietly = T))
suppressMessages(library("ggplot2", quietly = T))
suppressMessages(library("microbiome", quietly = T))
suppressMessages(library("DESeq2", quietly = T))
#######################################
alphabetaplot <- function(rare_data, bray_curtis, condition="condition"){
  # plot alpha and beta diversity of rarefied data
  p1 <- plot_richness(rare_data , measures = "Shannon", color = condition)+
    facet_wrap(condition, scales = "free_x", nrow = 1)+
    ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)
  p2 <- plot_richness(rare_data, x=condition, measures=c("Observed", "Shannon", "Chao1")) + geom_boxplot()
  p3 <- plot_ordination(rare_data, bray_curtis, color = condition)+stat_ellipse()+
    ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)
  
  suppressMessages(ggsave(paste0(name,"_alpha_all_samples.png"),p1))
  suppressMessages(ggsave(paste0(name,"_alpha_box_condition.png"),p2))
  suppressMessages(ggsave(paste0(name,"_beta_bray.png"),p3))
  
}
betaSplot <- function(rare_data, bray_curtis, condition="condition",p_name="p_name" ){
  p1 <- plot_ordination(rare_data, bray_curtis, color = p_name, shape = condition)+
    ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)+
    geom_line(aes(group=p_name))
  suppressMessages(ggsave(paste0(name,"_beta_sample_ordinate.png"),p1))
}
allgenuspercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Genus") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  suppressMessages(ggsave(paste0(name,"_Genus_barplot.png"),plot.mpg))
  suppressMessages(ggsave(paste0(name,"_Genus_barplot_legend.png"),legend,dpi = 300,width = 5000,height = 4000,units = "px"))
}
allorderpercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Order") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  suppressMessages(ggsave(paste0(name,"_Order_barplot.png"),plot.mpg))
  suppressMessages(ggsave(paste0(name,"_Order_barplot_legend.png"),legend,dpi = 300,width = 5000,height = 4000,units = "px"))
}
allfamilypercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Family") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  suppressMessages(ggsave(paste0(name,"_Family_barplot.png"),plot.mpg))
  suppressMessages(ggsave(paste0(name,"_Family_barplot_legend.png"),legend,dpi = 300,width = 5000,height = 4000,units = "px"))
}
allphylumpercondition <- function(rare_data, condition="condition"){
  bar_one <- plot_bar(rare_data, x=condition, fill="Phylum") + facet_wrap(~condition, scales="free_x") + 
    geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='bottom')
  suppressMessages(ggsave(paste0(name,"_Phylum_barplot.png"),plot.mpg))
}
bar_compare <- function(ps.top20, condition="condition", comp_name){
  p1 <- plot_bar(ps.top20, "condition", fill="condition", facet_grid=~Genus)+ 
    geom_bar(aes(color=condition, fill=condition), stat="identity", position="stack") + theme(legend.position="bottom")
  p2 <- plot_bar(ps.top20, "Genus", "Abundance", "Genus", 
           facet_grid="condition") + 
    geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
  suppressMessages(ggsave(paste0(name,"_bar_",comp_name, ".png"),p1))
  suppressMessages(ggsave(paste0(name,"_bar_",comp_name, "2.png"),p2))
}
DAVOLC <- function(pseq2, condition="condition", cond1="21-days", cond2="Before_treatment"){
  diagdds = phyloseq_to_deseq2(pseq2, ~condition)
  
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds, fitType="local")
  
  
  # investigate the results 
  res = results(diagdds,contrast = c(condition,cond1,cond2))
  res = res[order(res$padj, na.last=NA), ]
  alpha = 0.01
  sigtab = res[(res$padj < alpha), ]
  suppressMessages(library("EnhancedVolcano", quietly = T))
  svg(paste(name, "_volcanoplot_E_non_Filtered.svg", sep=""),  width=15, height=10)
  volc <- EnhancedVolcano(res,
                          title = "",
                          subtitle = bquote(italic(Volcanoplot)),
                          lab = rownames(res),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          xlab = bquote(~Log[2]~ 'fold change'),
                          pCutoff = 10e-4,
                          FCcutoff = 5,
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
############# for CL ##################
ps <- readRDS(opt$input)
rare_data <- rarefy_even_depth(ps, rngseed =1, sample.size = opt$rare)
pseq2 <- aggregate_taxa(rare_data, "Genus") 
pseq.rel <- microbiome::transform(rare_data, "compositional") # raltive counts 
# bar plot of top species per condition
top20 <- names(sort(taxa_sums(pseq2), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(pseq2, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
bray_curtis <- ordinate(rare_data, method = "PCoA")
bray_curtis_rel <- ordinate(pseq.rel, method = "PCoA")
target <- c("Lactobacillus","Streptococcus","Staphylococcus","Escherichia-Shigella",
            "Klebsiella")
ex3 <- subset_taxa(microbiome::transform(rare_data, "compositional"), Genus %in% target)
###################################
name=opt$name
setwd(opt$working_dir)
write(paste0("\033[0;", 32, "m","Alpha and Beta started...","\033[0m"), stderr())
alphabetaplot(rare_data, bray_curtis)
betaSplot(rare_data,bray_curtis)
write(paste0("\033[0;", 32, "m","Bar plots started...","\033[0m"), stderr())
allgenuspercondition(pseq.rel,condition = opt$condition)
allorderpercondition(pseq.rel,condition = opt$condition)
allfamilypercondition(pseq.rel,condition = opt$condition)
allphylumpercondition(pseq.rel,condition = opt$condition)
bar_compare(ps.top20,opt$condition,"top20")
bar_compare(ex3,opt$condition,"target_genus")
write(paste0("\033[0;", 32, "m","DA started...","\033[0m"), stderr())
DAVOLC(pseq2, condition = opt$condition, cond1 = opt$sample, cond2 = opt$control)