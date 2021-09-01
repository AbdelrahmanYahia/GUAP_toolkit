library("tidyverse")
library("qiime2R")
library("phyloseq")
library("ggplot2")
library("microbiome")

#############################################
# import data 
SVs<-read_qza("QIIME2/table.qza")
metadata<-read_q2metadata("sample-metadata.tsv")
taxonomy<-read_qza("QIIME2/classify/taxonomy.qza")

ps <- qza_to_phyloseq(
  features="QIIME2/table.qza",
  tree="QIIME2/align/rooted-tree.qza", "QIIME2/classify/taxonomy.qza",
  metadata = "sample-metadata.tsv"
)



psf <- subset_taxa(ps, !is.na(Genus) & !Genus %in% c("", "uncharacterized", NA, "unknow", "Unknown"))
pseq2 <- aggregate_taxa(psf, "Genus") 

rare_data <- rarefy_even_depth(ps, rngseed =1 , sample.size = 12000)

# plot alpha and beta diversity of rarefied data
plot_richness(rare_data , measures = "Shannon", color = "condition")+
  facet_wrap("condition", scales = "free_x", nrow = 1)+
  ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)

plot_richness(rare_data, x="condition", measures=c("Observed", "Shannon")) + geom_boxplot()

bray_curtis <- ordinate(rare_data, method = "PCoA")
jaccard <- ordinate(rare_data, method ="PCoA", distance = "jaccard")

plot_ordination(rare_data, bray_curtis, color = "condition")+stat_ellipse()+
  ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)
plot_ordination(rare_data, bray_curtis, color ="condition")+
  stat_ellipse()

p2 = plot_ordination(rare_data, bray_curtis, type="condition", color="condition") 
p2 + geom_polygon(aes(fill=condition)) + geom_point(size=2, alpha = 0.5) + ggtitle("samples")


# plot alpha and beta diversity of RAW data 
plot_richness(psf , measures = "Shannon", color = "condition",shape = "Gender")+
  facet_wrap("condition", scales = "free_x", nrow = 1)+
  ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)

plot_richness(psf, x="condition", measures=c("Observed", "Shannon")) + geom_boxplot()

bray_curtis <- ordinate(psf, method = "PCoA")
jaccard <- ordinate(psf, method ="PCoA", distance = "jaccard")

plot_ordination(psf, bray_curtis, color = "condition", shape = "Gender")+stat_ellipse()+
  ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)
plot_ordination(psf, bray_curtis, color = "condition")+
  stat_ellipse()

p2 = plot_ordination(psf, bray_curtis, type="condition", color="condition", shape="Gender") 
p2 + geom_polygon(aes(fill=condition)) + geom_point(size=2, alpha = 0.5) + ggtitle("samples")

# top 20 genus Heatmap 
top20g <- names(sort(taxa_sums(pseq2), decreasing=TRUE))[1:20]
ps.top20g <- transform_sample_counts(pseq2, function(OTU) OTU/sum(OTU))
ps.top20g <- prune_taxa(top20g, ps.top20g)
tablg <- (as.data.frame(otu_table(ps.top20g)))

p1 <-ComplexHeatmap::pheatmap(scale(tablg), annotation_col = metadata,
                              col = colorRampPalette(c("deepskyblue4", "mintcream"))(50),
                              column_split = metadata$condition,
                              cluster_rows = F, cluster_cols = F,
                              border_color = NA, scale = T)

# bar plot of top 20 taxa 
par(mar = c(13, 4, 4, 2) + 0.1) # make more room on bottom margin
N <- 20
barplot(sort(taxa_sums(pseq2), T)[1:N]/nsamples(pseq2), las=2)


metadata<-read_q2metadata("sample-metadata.tsv")
SVs<-read_qza("table.qza")$data
taxonomy<-read_qza("taxonomy.qza")$data

SVs<-apply(SVs, 2, function(x) x/sum(x)*100) #convert to percent

SVsToPlot<-  
  data.frame(MeanAbundance=rowMeans(SVs)) %>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(30, MeanAbundance) %>%
  pull(Feature.ID) #extract only the names from the table

SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  mutate(Feature.ID=if_else(Feature.ID %in% SVsToPlot,  Feature.ID, "Remainder")) %>% #flag features to be collapsed
  group_by(SampleID, Feature.ID) %>%
  summarize(Abundance=sum(Abundance)) %>%
  left_join(metadata) %>%
  mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature=paste(Feature.ID, Taxon)) %>%
  mutate(Feature=gsub("[kpcofgs]__", "", Feature)) %>% # trim out leading text from taxonomy string
  ggplot(aes(x=SampleID, y=Feature, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`Sample-Sample-condition`, scales="free_x") +
  theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)")
results<-read_qza("differentials.qza")$data

results<-results %>% mutate(Significant=if_else(we.eBH<0.1,"*", ""))

tree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% results$Feature.ID]) # remove all the features from the tree we do not have data for
ggtree(tree, layout="circular") %<+% results +
  geom_tippoint(aes(fill=diff.btw), shape=21, color="grey50")  +
  geom_tiplab2(aes(label=Significant), size=10) +
  scale_fill_gradient2(low="darkblue",high="darkred", midpoint = 0, mid="white", name="log2(fold-change") +
  theme(legend.position="right")











library(tidyverse)
library(qiime2R)

metadata<-read_q2metadata("sample-metadata.tsv")
uwunifrac<-read_qza("unweighted_unifrac_pcoa_results.qza")
shannon<-read_qza("shannon_vector.qza")$data %>% rownames_to_column("SampleID") 

uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes(x=PC1, y=PC2, color=`body-site`, shape=`reported-antibiotic-usage`, size=shannon)) +
  geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_shape_manual(values=c(16,1), name="Antibiotic Usage") + #see http://www.sthda.com/sthda/RDoc/figure/graphs/r-plot-pch-symbols-points-in-r.png for numeric shape codes
  scale_size_continuous(name="Shannon Diversity") +
  scale_color_discrete(name="Body Site")
ggsave("PCoA.pdf", height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches

plot_tree(ps, color = "condition", label.tips = "Phylum", size = "abundance", plot.margin = 0.5, ladderize = TRUE,text.size = 3)


prp11 <- ggplot(df, aes(fill=OTU, y=Abundance, x=Sample)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  theme_ipsum() +
  xlab("")
legend <- cowplot::get_legend(prp11)
plot.mpg <- prp11 + theme(legend.position='none')
print(plot.mpg)
# ggsave(prp11, file='prp11.pdf',device = cairo_pdf)

