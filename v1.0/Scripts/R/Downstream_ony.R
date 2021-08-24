library("dada2")
library("phyloseq")
library("Biostrings")
library("ggplot2")
library("vegan")
library("decontam")
library("knitr")
library("BiocStyle")
library("DECIPHER")
library("phangorn")
library("gridExtra")
library("viridis")
library("hrbrthemes")
library("tidyverse")
library("metagMisc")
library("ComplexHeatmap")
library("microbiome")
library("eulerr")
library("microbiomeutilities")

load("AUC_all.RData")

################       downstream analysis    ##################################
###### check counts #####
t.track <- t(track)
t.track <- as.data.frame(t.track)
print((t.track) %>% 
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
        ylab(""))
t.track2 <- as.data.frame(t(data.frame(
  input=track$input,
  filttered=track$nonchim
)))

names(t.track2) <- names(t.track)
p11 <- (t.track2) %>% 
  pivot_longer(cols = names((t.track)), 
               names_to = "gene", 
               values_to = "expression")

p11$data <- c(rep(c("input"), times = 18),rep(c("output"), times = 18))

ggplot(p11, aes(fill=data, y=expression, x=gene)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Counts filttered") +
  theme_ipsum() +
  xlab("")
################################################################################
# load sample sheet 
samdf <- read.csv("samples.csv")
rownames(samdf) <- samdf$sample.name

# create phyloseq obj
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# filter NA till genus level
psf <- subset_taxa(ps, !is.na(Genus) & !Genus %in% c("", "uncharacterized", NA, "unknow", "Unknown"))

# aggregate taxa @ genus level
pseq2 <- aggregate_taxa(psf, "Genus") 

# rarefication
rarecurve(t(otu_table(ps)),label = F)
rare_data <- rarefy_even_depth(ps, rngseed =1 )

# plot alpha and beta diversity of rarefied data
plot_richness(rare_data , measures = "Shannon", color = "condition",shape = "Gender")+
  facet_wrap("condition", scales = "free_x", nrow = 1)+
  ggtitle("Alpha diversity (Shannon)")+ geom_point(size=2, alpha=0.7)

plot_richness(rare_data, x="condition", measures=c("Observed", "Shannon")) + geom_boxplot()

bray_curtis <- ordinate(rare_data, method = "PCoA")
jaccard <- ordinate(rare_data, method ="PCoA", distance = "jaccard")

plot_ordination(rare_data, bray_curtis, color = "condition", shape = "Gender")+stat_ellipse()+
  ggtitle("Beta diversity (Bray-Curtis PCoA)")+ geom_point(size=3, alpha=1)
plot_ordination(rare_data, bray_curtis, color = "condition")+
  stat_ellipse()

p2 = plot_ordination(rare_data, bray_curtis, type="condition", color="condition", shape="Gender") 
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

p1 <-ComplexHeatmap::pheatmap(scale(tablg), annotation_col = samdf,
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
                         prevalence = 0.75, include.lowest=F)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}
mycols <- c(nonCRC="#d6e2e9", CRC="#cbf3f0", H="#fcf5c7") 
plot(venn(list_core),
     fills = mycols)
taxa_names(pseq.rel)[1:10]


# bar plot of normalized per condition 
bar_one <- plot_bar(microbiome::transform(psf, "compositional"), x="condition", fill="Genus") + facet_wrap(~condition, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") + theme(legend.position="bottom")
legend <- cowplot::get_legend(bar_one)
plot.mpg <- bar_one + theme(legend.position='none')
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
plot_landscape(psf, "NMDS", "bray", col = "condition")
