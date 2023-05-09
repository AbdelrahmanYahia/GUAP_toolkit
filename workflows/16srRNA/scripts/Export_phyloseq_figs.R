library("optparse", quietly = T)
defaultW <- getOption("warn") 
options(warn = -1) 
#######################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="phyloseq object RDS", metavar="character"),
  make_option(c("-n", "--name"), type="character", default="GUAP-DE-out", 
              help="Name of analysis [default= %default]", metavar="character"),
  make_option(c("-c", "--condition"), type="character", default="condition", 
              help="Condition column name to generate figures [default= %default]", metavar="character"),
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
suppressMessages(library("data.table", quietly = T))
suppressMessages(library("DESeq2", quietly = T))
suppressMessages(library("ggpubr", quietly = T))
suppressMessages(library("microbiomeutilities", quietly = T))
suppressMessages(library("picante", quietly = T))
suppressMessages(library("vegan", quietly = T))
suppressMessages(library("RColorBrewer", quietly = T))
#######################################
alphabetaplot <- function(rare_data, bray_curtis, condition){
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
betaSplot <- function(rare_data, bray_curtis, condition,p_name ){
  p1 <- plot_ordination(rare_data, bray_curtis, color = p_name, shape = condition)+
    ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)+
    geom_line(aes(group=p_name))
  suppressMessages(ggsave(paste0(name,"_beta_sample_ordinate.png"),p1))
}
allgenuspercondition <- function(rare_data, condition){
  bar_one <- plot_bar(rare_data, x=condition, fill="Genus") + facet_wrap(facets = c(condition), scales="free_x") + 
    geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  suppressMessages(ggsave(paste0(name,"_Genus_barplot.png"),plot.mpg))
  suppressMessages(ggsave(paste0(name,"_Genus_barplot_legend.png"),legend,dpi = 300,width = 6000,height = 5000,units = "px"))
}
allorderpercondition <- function(rare_data, condition){
  bar_one <- plot_bar(rare_data, x=condition, fill="Order") + facet_wrap(facets = c(condition), scales="free_x") +  # nolint
    geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  suppressMessages(ggsave(paste0(name,"_Order_barplot.png"),plot.mpg))
  suppressMessages(ggsave(paste0(name,"_Order_barplot_legend.png"),legend,dpi = 300,width = 5000,height = 4000,units = "px"))
}
allfamilypercondition <- function(rare_data, condition){
  bar_one <- plot_bar(rare_data, x=condition, fill="Family") + facet_wrap(facets = c(condition), scales="free_x") + 
    geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='none')
  suppressMessages(ggsave(paste0(name,"_Family_barplot.png"),plot.mpg))
  suppressMessages(ggsave(paste0(name,"_Family_barplot_legend.png"),legend,dpi = 300,width = 5000,height = 4000,units = "px"))
}
allphylumpercondition <- function(rare_data, condition){
  bar_one <- plot_bar(rare_data, x=condition, fill="Phylum") + facet_wrap(facets = c(condition), scales="free_x") + 
    geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme(legend.position="bottom")
  legend <- cowplot::get_legend(bar_one)
  plot.mpg <- bar_one + theme(legend.position='bottom')
  suppressMessages(ggsave(paste0(name,"_Phylum_barplot.png"),plot.mpg,dpi = 300,width = 5000,height = 4000,units = "px"))
}
bar_compare <- function(ps.top20, condition, comp_name){
  p1 <- plot_bar(ps.top20, condition, fill=condition, facet_grid=~Genus)+ 
    geom_bar(aes(color=condition, fill=condition), stat="identity", position="stack") + theme(legend.position="bottom")
  p2 <- plot_bar(ps.top20, "Genus", "Abundance", "Genus", 
           facet_grid=condition) + 
    geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
  suppressMessages(ggsave(paste0(name,"_bar_",comp_name, ".png"),p1, dpi = 300,width = 5000,height = 4000,units = "px"))
  suppressMessages(ggsave(paste0(name,"_bar_",comp_name, "2.png"),p2, dpi = 300,width = 5000,height = 4000,units = "px"))
}
plotheatmap <- function(ps_new, level, n){

  ps_new_agg <- aggregate_taxa(ps_new,level)
  # heatmap 
  top20g <- names(sort(taxa_sums(ps_new_agg), decreasing=TRUE))[1:30]
  ps.top20d <- transform_sample_counts(ps_new_agg, function(OTU) OTU/sum(OTU))
  ps.top20d <- prune_taxa(top20g, ps.top20d)
  hm1 <- plot_heatmap(ps.top20d)

  top20g <- names(sort(taxa_sums(ps_new_agg), decreasing=TRUE))[1:250]
  ps.top20d <- transform_sample_counts(ps_new_agg, function(OTU) OTU/sum(OTU))
  ps.top20d <- prune_taxa(top20g, ps.top20d)
  hm2 <- plot_heatmap(ps.top20d)
  suppressMessages(pdf(paste(name, "_",level,"_Heatmap1.pdf", sep=""),  width=15, height=20))
  print(hm1)
  suppressMessages(suppressMessages(dev.off()))
  suppressMessages(pdf(paste(name, "_",level,"_Heatmap2.pdf", sep=""),  width=15, height=25))
  print(hm2)
  suppressMessages(suppressMessages(dev.off()))
}

plots_with_sign <- function (ps, cond){
  print_rez <- function(res,condition="condition",cond1="2_7-days",cond2="3_21-days",pcutoff = 0.01,genre ){
    res = results(diagdds, cooksCutoff = FALSE,contrast = c(condition,cond2,cond1))
    res.df <- as.data.frame(res)
    write.csv(res.df, paste0(cond1,"vs",cond2,"_",genre,"_res.csv"), quote = F)
    res.df$X <- rownames(res.df)
    res.df = res.df[order(res.df$padj, na.last=NA), ]
    alpha = pcutoff
    sigtab = res.df[(res.df$padj < alpha), ]
    sigtab$X = factor(sigtab$X,levels=sigtab$X[order(sigtab$log2FoldChange)])
    
    p = ggplot(data = sigtab, 
              aes(x = X, y =  log2FoldChange, fill =  log2FoldChange, color =  log2FoldChange)) + 
      geom_bar(stat = "identity", width = 0.7, 
              position = position_dodge(width = 0.4)) +
      geom_errorbar(aes(ymin =  log2FoldChange -  lfcSE, ymax = log2FoldChange +  lfcSE), width = 0.2,
                    position = position_dodge(0.05), color = "black") + 
      labs(x = NULL, y = "Log fold change", 
          title = paste(cond1,"vs",cond2)) + 
      theme_bw() + 
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            panel.grid.minor.y = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1))
    return(p)
  }
  plot_bta_box_condition <- function(ps1.rel,wt.unifrac.dist.w){
    sub_dist <- list()
    groups_all <- sample_data(ps1.rel)$condition
    abrel_bray <- as.matrix(wt.unifrac.dist.w)
    for (group in levels(groups_all)) { 
      row_group <- which(groups_all == group)
      sample_group <- sample_names(ps1.rel)[row_group]
      sub_dist[[group]] <- abrel_bray[ sample_group, sample_group]
      sub_dist[[group]][!lower.tri(sub_dist[[group]])] <- NA
    }
    braygroups<- melt(sub_dist)
    df.bray <- braygroups[complete.cases(braygroups), ]
    df.bray$L1 <- factor(df.bray$L1, levels=names(sub_dist))
    
    conditions <- levels(ps1.metadata$condition)
    conditions.pairs <- combn(seq_along(conditions), 2, simplify = FALSE, FUN = function(i)conditions[i])
    
    ppppp <- ggboxplot(df.bray, x = "L1", y = "value")+
      stat_compare_means(comparisons = conditions.pairs)+ # Add significance levels
      stat_compare_means(label.y = max(df.bray$value)*1.4) 
    return(ppppp)
  }
  ps1 <- prune_taxa(taxa_sums(ps) > 1, ps) 
  metadata <- sample_data(ps1)
  new_metadata <- data.frame(
    IDs = rownames(metadata),
    condition = metadata[,cond]
  )
  names(new_metadata) <- c("IDs", "condition")
  sample_data(ps1) <- new_metadata

  ##############################################################
  ###############             Stats          ##################
  ##############################################################

  reads_per_OTU <- taxa_sums(ps1)
  ps1.dt.taxa = data.table(tax_table(ps1),OTUabundance = taxa_sums(ps1),
                          OTU = taxa_names(ps1))
  n_reads <- sum(reads_per_OTU)
  n_ASV <- ntaxa(ps1)
  n_10_ASV <- length(reads_per_OTU[reads_per_OTU < 10])
  n_10_R <- sum(reads_per_OTU[reads_per_OTU < 10])
  per_10_R <- round(((sum(reads_per_OTU[reads_per_OTU < 10])/sum(reads_per_OTU))*100),2)
  per_10_ASV <- round((n_10_ASV/n_ASV*100),2)
  n_singltones <- ps1.dt.taxa[(OTUabundance <= 0), .N]
  n_doubles <- ps1.dt.taxa[(OTUabundance <= 2), .N]
  per_n_double_ASV <- round((n_doubles/n_ASV*100),3)
  per_n_double_R <- round((sum(reads_per_OTU[reads_per_OTU <= 2])/sum(reads_per_OTU)*100),3)

  ps_info <- data.frame(
    "Number of reads" = n_reads,
    "Number of ASVs" = n_ASV,
    "ASVs with less than 10 reads" = n_10_ASV,
    "ASVs with less than 10 reads %" = per_10_ASV,
    "total number of reads in ASVs with less than 10 reads" = n_10_R,
    "total number of reads in ASVs with less than 10 reads %" = per_10_R,
    "Number of singletones" = n_singltones,
    "Number of Doubles (ASVs)" = n_doubles,
    "Number of Doubles (ASVs) %" = per_n_double_ASV,
    "Number of Doubles (Reads %)" = per_n_double_R
  )
  write(paste0("\033[0;", 32, "m","exporting ps info...","\033[0m"), stderr())

  write.csv(ps_info, paste0("ps_info.csv"), quote = F)



  ##############################################################
  ###############       Alpha Diversity       ##################
  ##############################################################

  otu_table_ps1 <- as.data.frame(ps1@otu_table)
  metadata_table_ps1  <- as.data.frame(ps1@sam_data)
  df.pd <- pd(t(otu_table_ps1), ps1@phy_tree,include.root=F)
  metadata_table_ps1$Phyogenetic_diversity <- df.pd$PD 
  ps1.adiv <- estimate_richness(ps1, measures = c("Chao1", 
                                                  "Shannon", 
                                                  "Observed", 
                                                  "InvSimpson"))
  ps1.metadata <- as(sample_data(ps1), "data.frame")
  ps1.metadata$Observed <- ps1.adiv$Observed 
  ps1.metadata$Shannon <- ps1.adiv$Shannon
  ps1.metadata$InvSimpson <- ps1.adiv$InvSimpson
  ps1.metadata$Phyogenetic_diversity <- df.pd$PD
  conditions <- levels(ps1.metadata$condition)
  conditions.pairs <- combn(seq_along(conditions), 2, simplify = FALSE, FUN = function(i)conditions[i])

  # observed 
  ov <- ggviolin(ps1.metadata, x = "condition", y = "Observed", fill = "condition",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = conditions.pairs)+ # Add significance levels
    stat_compare_means(label.y = max(ps1.metadata$Observed)*1.4)                                      

  # Shannon 
  sv <- ggviolin(ps1.metadata, x = "condition", y = "Shannon", fill = "condition",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = conditions.pairs)+ # Add significance levels
    stat_compare_means(label.y = max(ps1.metadata$Shannon)*1.4)  

  # Phyogenetic_diversity 
  pv <- ggviolin(ps1.metadata, x = "condition", y = "Phyogenetic_diversity", fill = "condition",
          palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = conditions.pairs)+ # Add significance levels
    stat_compare_means(label.y = max(ps1.metadata$Phyogenetic_diversity)*1.4)                                      
  write(paste0("\033[0;", 32, "m","Ploting alpha diversity...","\033[0m"), stderr())
  suppressMessages(ggsave(paste0(name,"_observed.png"),ov))
  suppressMessages(ggsave(paste0(name,"_shannon.png"),sv))
  suppressMessages(ggsave(paste0(name,"_Phyogenetic_diversity.png"),pv))
  
  ##############################################################
  ###############        beta Diversity       ##################
  ##############################################################
  write(paste0("\033[0;", 32, "m","plotting beta diversity...","\033[0m"), stderr())
  set.seed(19773)

  ps1.rel <- microbiome::transform(ps1, "relative.abundance")
  metadf <- data.frame(sample_data(ps1.rel))

  # weighted unifrac relative abundance
  wt.unifrac.dist.w <- phyloseq::distance(ps1.rel,method = "wunifrac")
  adonis.wt.unifrac.dist.w <- vegan::adonis2(wt.unifrac.dist.w~ condition, metadf)

  wt.unifrac.dispr <- vegan::betadisper(wt.unifrac.dist.w, 
                                        phyloseq::sample_data(ps1.rel)$condition)
  permanova.wt.unifrac.dist.w <- permutest(wt.unifrac.dispr)
  write.csv(permanova.wt.unifrac.dist.w$tab, paste0("permanova.wt.unifrac.dist.w.csv"), quote = F)
  write.csv(adonis.wt.unifrac.dist.w, paste0("adonis.wt.unifrac.dist.w.csv"), quote = F)

  ordu.wt.uni = ordinate(ps1.rel, "PCoA", "wunifrac")
  wt.unifrac <- plot_ordination(ps1.rel, ordu.wt.uni, color = "condition") 
  wt.unifrac <- wt.unifrac + 
    ggtitle("Weighted UniFrac") + 
    geom_point(size = 3)
  suppressMessages(png(filename=paste(name, "_","PCoA_wt_unifrac.png", sep=""), res=150, width=7, height= 5, units= "in"))
  print(wt.unifrac)
  suppressMessages(dev.off())
  suppressMessages(png(filename=paste(name, "_","box_plot_wt_unifrac.png", sep=""), res=150, width=7, height= 5, units= "in"))
  boxplot(wt.unifrac.dispr, main = "", xlab = "")
  suppressMessages(dev.off())
  suppressMessages(png(filename=paste(name, "_","box_with_sign_wt_unifrac.png", sep=""), res=150, width=7, height= 5, units= "in"))
  print(plot_bta_box_condition(ps1.rel,wt.unifrac.dist.w))
  suppressMessages(dev.off())
  # unifrac relative abundance
  unifrac.dist <- phyloseq::distance(ps1.rel,method = "unifrac")
  adonis.unifrac.dist <- vegan::adonis2(unifrac.dist~ condition, metadf)
  unifrac.dispr <- vegan::betadisper(unifrac.dist, 
                                    phyloseq::sample_data(ps1.rel)$condition)
  permanova.unifrac.dist <- permutest(unifrac.dispr)
  write.csv(permanova.unifrac.dist$tab, paste0("permanova.unifrac.dist.csv"), quote = F)
  write.csv(adonis.unifrac.dist, paste0("adonis.unifrac.dist.csv"), quote = F)
  ordu.uni = ordinate(ps1.rel, "PCoA", "unifrac")
  unifrac <- plot_ordination(ps1.rel, ordu.uni, color = "condition") 
  unifrac <- unifrac + 
    ggtitle("UniFrac") + 
    geom_point(size = 3)
  suppressMessages(png(filename=paste(name, "_","PCoA_unifrac.png", sep=""), res=150, width=7, height= 5, units= "in"))
  print(unifrac)
  suppressMessages(dev.off())
  suppressMessages(png(filename=paste(name, "_","box_plot_unifrac.png", sep=""), res=150, width=7, height= 5, units= "in"))
  boxplot(unifrac.dispr, main = "", xlab = "")
  suppressMessages(dev.off())
  suppressMessages(png(filename=paste(name, "_","box_with_sign_unifrac.png", sep=""), res=150, width=7, height= 5, units= "in"))
  print(plot_bta_box_condition(ps1.rel,unifrac.dist))
  suppressMessages(dev.off())
  # bray relative abundance
  bray.dist <- phyloseq::distance(ps1.rel,method = "bray")
  adonis.bray.dist <- vegan::adonis2(bray.dist~ condition, metadf)
  bray.dispr <- vegan::betadisper(bray.dist, 
                                    phyloseq::sample_data(ps1.rel)$condition)
  permanova.bray.dist <- permutest(bray.dispr)
  write.csv(permanova.bray.dist$tab, paste0("permanova.bray.dist.csv"), quote = F)
  write.csv(adonis.bray.dist, paste0("adonis.bray.dist.csv"), quote = F)
  ordu.bray = ordinate(ps1.rel, "PCoA", "bray")
  bray <- plot_ordination(ps1.rel, ordu.bray, color = "condition") 
  bray <- bray + 
    ggtitle("bray") + 
    geom_point(size = 3)

  suppressMessages(png(filename=paste(name, "_","PCoA_bray.png", sep=""), res=150, width=7, height= 5, units= "in"))
  print(bray)
  suppressMessages(dev.off())
  suppressMessages(png(filename=paste(name, "_","box_plot_bray.png", sep=""), res=150, width=7, height= 5, units= "in"))
  boxplot(bray.dispr, main = "", xlab = "")
  suppressMessages(dev.off())
  suppressMessages(png(filename=paste(name, "_","box_with_sign_bray.png", sep=""), res=150, width=7, height= 5, units= "in"))
  print(plot_bta_box_condition(ps1.rel,bray.dist))
  suppressMessages(dev.off())

  ### test with geomtric means ###
  write(paste0("\033[0;", 32, "m","DA Genus...","\033[0m"), stderr())
  pseq2 <- aggregate_taxa(ps1, "Genus") 

  diagdds = phyloseq_to_deseq2(pseq2, ~condition)
  gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds)
  res = results(diagdds)

  conditions <- levels(ps1.metadata$condition)
  conditions.pairs <- combn(seq_along(conditions), 2, simplify = FALSE, FUN = function(i)conditions[i])
  for ( i in conditions.pairs ){
    one <- (i[1])
    two <- (i[2])
    suppressMessages(png(filename=paste(two, "VS",one,"_","DA_Genus.png", sep=""), res=150, width=7, height= 5, units= "in"))
    print(print_rez(res, "condition", one, two,genre="Genus"))
    suppressMessages(dev.off())
  }

  write(paste0("\033[0;", 32, "m","DA Phylum...","\033[0m"), stderr())
  pseq2 <- aggregate_taxa(ps1, "Phylum") 

  diagdds = phyloseq_to_deseq2(pseq2, ~condition)
  gm_mean = function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  diagdds = DESeq(diagdds)
  res = results(diagdds)

  conditions <- levels(ps1.metadata$condition)
  conditions.pairs <- combn(seq_along(conditions), 2, simplify = FALSE, FUN = function(i)conditions[i])
  for ( i in conditions.pairs ){
    one <- (i[1])
    two <- (i[2])
    suppressMessages(png(filename=paste(two, "VS",one,"_","DA_Phylum.png", sep=""), res=150, width=7, height= 5, units= "in"))
    print(print_rez(res, "condition", one, two,genre="Phylum"))
    suppressMessages(dev.off())
  }

}


ps <- readRDS(opt$input)
min_sample_size <- min(sample_sums(ps))

suppressMessages(rare_data <- rarefy_even_depth(ps, rngseed =1, sample.size = min_sample_size))
pseq2 <- aggregate_taxa(rare_data, "Genus") 
pseq.rel <- microbiome::transform(rare_data, "compositional") # raltive counts 
# bar plot of top species per condition
top20 <- names(sort(taxa_sums(pseq2), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(pseq2, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
bray_curtis <- ordinate(rare_data, method = "PCoA")
bray_curtis_rel <- ordinate(pseq.rel, method = "PCoA")

###################################

name=opt$name
setwd(opt$working_dir)

tryCatch({
    write(paste0("\033[0;", 32, "m","Plots started...","\033[0m"), stderr())
    plots_with_sign(ps, opt$condition)
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error in Plots with significance: ","\033[0m",err), stderr())
  }
)

write(paste0("\033[0;", 32, "m","Alpha and Beta started...","\033[0m"), stderr())

tryCatch({
    richness_RAW <- (estimate_richness(ps, split = TRUE, measures = NULL))
    write.csv(richness_RAW, paste0("Richness_RAW.csv"), quote = F)
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error in estimate RAW richness: ","\033[0m",err), stderr())
  }
)
tryCatch({
    richness_relative <- (estimate_richness(rare_data, split = TRUE, measures = NULL))
    write.csv(richness_relative, paste0("richness_rarified.csv"), quote = F)
  }, error = function(err){
    write(paste0("\033[0;", 31, "m","Error in estimate rarified richness: ","\033[0m",err), stderr())
  }, finally = {
    alphabetaplot(rare_data, bray_curtis, opt$condition)
  }
)

tryCatch({
    write(paste0("\033[0;", 32, "m","Bar plots started...","\033[0m"), stderr())
    allgenuspercondition(pseq.rel,condition = opt$condition)
    allorderpercondition(pseq.rel,condition = opt$condition)
    allfamilypercondition(pseq.rel,condition = opt$condition)
    allphylumpercondition(pseq.rel,condition = opt$condition)
    bar_compare(ps.top20,opt$condition,"top20")
  }, error = function(err){
      write(paste0("\033[0;", 31, "m","Error in Plotting bars: ","\033[0m",err), stderr())
  }
)

tryCatch({
    write(paste0("\033[0;", 32, "m","Exporting Heatmaps...","\033[0m"), stderr())
    plotheatmap(ps, "Genus", 30)
    plotheatmap(ps, "Order", 30)
    plotheatmap(ps, "Family", 30)
    plotheatmap(ps, "Species", 30)
  }, error = function(err){
      write(paste0("\033[0;", 31, "m","Error in Plotting heatmaps: ","\033[0m",err), stderr())
  }
)

tryCatch({
    write(paste0("\033[0;", 32, "m","Exporting RAW Taxonomy tables...","\033[0m"), stderr())
    dir.create("../Taxonomy_tables")
    setwd("../Taxonomy_tables")
    testnames <- as.data.frame(ps@tax_table@.Data)
    for ( i in names(testnames[1:7])){
      write.csv(as.data.frame(otu_table(aggregate_taxa(ps, i))), paste0(i,"_taxonomy_table.csv"), quote = F)
    }

    write(paste0("\033[0;", 32, "m","Exporting Relative Abundance (Rarefied) Taxonomy tables...","\033[0m"), stderr())
    dir.create("./Relative_Taxonomy_tables")
    setwd("./Relative_Taxonomy_tables")
    testnames <- as.data.frame(ps@tax_table@.Data)
    ps2 <- microbiome::transform(ps, "compositional") # raltive counts
    for ( i in names(testnames[1:7])){
      write.csv(as.data.frame(otu_table(aggregate_taxa(ps2, i))), paste0(i,"_taxonomy_table.csv"), quote = F)
    }
  }, error = function(err){
      write(paste0("\033[0;", 31, "m","Error in exporting abundace tables: ","\033[0m",err), stderr())
  }
)

write(paste0("\033[0;", 32, "m","Figures generated...","\033[0m"), stderr())