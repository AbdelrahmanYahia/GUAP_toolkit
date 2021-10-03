library("dada2")
library("phyloseq")
library("Biostrings")
library("ggplot2")
library("vegan")
library("knitr")
library("DECIPHER")
library("phangorn")
library("gridExtra")
library("viridis")
library("hrbrthemes")
library("tidyverse")
library("ComplexHeatmap")
library("microbiome")
library("eulerr")
library("microbiomeutilities")
library("qiime2R")
library("ggtree")
library("ggtreeExtra")

#############################################
# import data ls 
SVs<-read_qza("QIIME2/table.qza")

taxonomy<-read_qza("QIIME2/classify/taxonomy.qza")
tree <- read_qza("QIIME2/align/rooted-tree.qza")$data

ps <- qza_to_phyloseq(
  features="QIIME2/table.qza",
  tree="QIIME2/align/rooted-tree.qza", "QIIME2/classify/taxonomy.qza",
  metadata = "sample-metadata.tsv"
)

ps <- phyloseq(ps@otu_table, ps@tax_table, sample_data(metadata), read_tree(tree))

############# from downstream only #############

# filter NA till genus level
# aggregate taxa @ genus level

# rarefication
# rarecurve(t(otu_table(ps)),label = F)
rare_data <- rarefy_even_depth(ps, rngseed =1, sample.size = 4000)
pseq2 <- aggregate_taxa(rare_data, "Genus") 

GP <- prune_taxa(taxa_sums(rare_data) > 25, rare_data)
mergedGP <- merge_samples(GP, "condition")
mergedGP <- tax_glom(mergedGP,"Genus")
mergedGP@sam_data$condition <- c("Control_Acute", "Control_Chronic", "STZ_Acute", "STZ_Chronic")
plot_tree(mergedGP, color="condition",base.spacing=0.006,nodelabf=nodeplotblank, label.tips="Genus",text.size = 2)



melt_simple <- psmelt(mergedGP) %>%
  filter(Abundance < 1000000) %>%
  select(OTU, val=Abundance)

p <- ggtree(mergedGP, layout="fan", open.angle=10) + 
  geom_tippoint(mapping=aes(color=Genus), 
                size=1.5,
                show.legend=FALSE)
p <- rotate_tree(p, -90)

p <- p +
  geom_fruit(
    data=melt_simple,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,
      group=label,
      fill=Genus,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 

p <- p +
  scale_fill_discrete(
    name="Genus",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), # The title of legend 
    legend.text=element_text(size=7) # The label text of legend, the sizes should be adjust with dpi.
  )
p



# plot alpha and beta diversity of rarefied data
plot_richness(rare_data , measures = "Shannon", color = "condition")+
  facet_wrap("condition", scales = "free_x", nrow = 1)+
  ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)

plot_richness(rare_data, x="condition", measures=c("Observed", "Shannon")) + geom_boxplot()

bray_curtis <- ordinate(rare_data, method = "PCoA")
jaccard <- ordinate(rare_data, method ="PCoA", distance = "jaccard")

plot_ordination(rare_data, bray_curtis, color = "condition")+stat_ellipse()+
  ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)

plot_ordination(rare_data, bray_curtis, color = "p_name", shape = "condition",type="split")+
  ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)+
  geom_line(aes(group=p_name))

plot_ordination(rare_data, bray_curtis, type="split", color="Phylum", 
                shape="condition", title="biplot", label = "station") +  
  geom_point(size=3)
# p2 = plot_ordination(rare_data, bray_curtis, color="condition") 
# p2 + geom_polygon(aes(fill=condition)) + geom_point(size=2, alpha = 0.5) + ggtitle("samples")


# plot alpha and beta diversity of RAW data 
plot_richness(ps , measures = "Shannon", color = "condition")+
  facet_wrap("condition", scales = "free_x", nrow = 1)+
  ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)

plot_richness(ps, x="condition", measures=c("Observed", "Shannon")) + geom_boxplot()

bray_curtis <- ordinate(ps, method = "PCoA")
jaccard <- ordinate(ps, method ="PCoA", distance = "jaccard")

plot_ordination(ps, bray_curtis, color = "condition")+stat_ellipse()+
  ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)


# p2 = plot_ordination(ps, bray_curtis, color="condition") 
# p2 + geom_polygon(aes(fill=condition)) + geom_point(size=2, alpha = 0.5) + ggtitle("samples")

# top 20 genus Heatmap 
top20g <- names(sort(taxa_sums(pseq2), decreasing=TRUE))[1:20]
ps.top20g <- transform_sample_counts(pseq2, function(OTU) OTU/sum(OTU))
ps.top20g <- prune_taxa(top20g, ps.top20g)
tablg <- (as.data.frame(otu_table(ps.top20g)))
samdf <- (ps.top20g@sam_data)
names(samdf) <-  c("sample_names","condition","pnames")
annotation <- as.data.frame(samdf$condition)
rownames(annotation) <- rownames(samdf)
p1 <-ComplexHeatmap::pheatmap(scale(tablg), annotation_col = annotation,
                              col = colorRampPalette(c("deepskyblue4", "mintcream"))(50),
                              column_split = samdf$condition,
                              cluster_rows = F, cluster_cols = F,
                              border_color = NA, scale = T)

# bar plot of top 20 taxa 
par(mar = c(13, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 20
barplot(sort(taxa_sums(pseq2), T)[1:N]/nsamples(pseq2), las=2)

# bar plot of raredata top 50 genus 
colnames(tax_table(rare_data)) <- c("Kingdom", "Phylum", "Class", "Order", "Family",  "Genus","Species")
genus <- tax_glom(rare_data, "Genus")
genus2 <- transform_sample_counts(genus,function(x)100*x/sum(x))
genus_sums <- names(sort(taxa_sums(genus2),TRUE)[1:50])
genus3 <- prune_taxa(genus_sums, genus2)
plot_bar(genus3, fill = "Genus")+geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")+
  facet_wrap("condition", scales = "free_x", nrow = 1)+
  ggtitle("Top 50 genera")

# venn diagram 

pseq.rel <- microbiome::transform(pseq2, "compositional") # raltive counts 
disease_states <- unique(as.character(meta(pseq.rel)$condition))
list_core <- c() # an empty object to store information

for (n in disease_states){ # for each variable n in DiseaseState
  print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, condition == n) # Choose sample from DiseaseState by n
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in at least 90% samples 
                         prevalence = 0.35, include.lowest=F)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)
taxa_names(pseq.rel)[1:10]
# 
# bef <- list_core$Before_treatment
# sev <- list_core$`7-days`
# twent <- list_core$`21-days`
# 
# Reduce(intersect, list(bef,sev,twent))
# 
# setdiff(bef, c(sev,twent))
# setdiff(sev, c(bef,twent))
# setdiff(twent, c(sev,bef))





# bar plot of normalized per condition 
bar_one <- plot_bar(microbiome::transform(rare_data, "compositional"), x="condition", fill="Genus") + facet_wrap(~condition, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
legend <- cowplot::get_legend(bar_one)
plot.mpg <- bar_one + theme(legend.position='none')
print(plot.mpg)
grid.newpage()
grid.draw(legend)

# bar plot of normalized per condition 
bar_one <- plot_bar(microbiome::transform(rare_data, "compositional"), x="condition", fill="Phylum") + facet_wrap(~condition, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") + theme(legend.position="bottom")
legend <- cowplot::get_legend(bar_one)
plot.mpg <- bar_one + theme(legend.position='bottom')
print(plot.mpg)
grid.newpage()
grid.draw(legend)


# bar plot of normalized samples
bar_one <- plot_bar(microbiome::transform(pseq2, "compositional"), fill="Genus") +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
legend <- cowplot::get_legend(bar_one)
plot.mpg <- bar_one + theme(legend.position='none')
print(plot.mpg)


# bar plot of top species per condition
top20 <- names(sort(taxa_sums(pseq2), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(pseq2, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

bar_one <- plot_bar(ps.top20, x="Genus",fill="Genus", facet_grid = ~condition) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")
legend <- cowplot::get_legend(bar_one)
plot.mpg <- bar_one + theme(legend.position='none')
print(plot.mpg)


# abundance 
plot_landscape(ps, "NMDS", "bray", col = "condition")

save.image("QIIME2_import_viz.RData")
################################################