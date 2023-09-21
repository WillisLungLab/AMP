library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(Maaslin2)

########################## Prepare Data ##########################

###Set Colors###
betadiv_colors <- c("#EF94A2", "#FB0207")
deseq2_colors <- c("#FB0207","#EF94A2")

###Load Data and Arrange###
#Load Oxygen OTU table
pankaj2 <- pankaj2 %>%
  tibble::column_to_rownames("Name")

#Load Oxygen Taxonomy Table
pankaj2.taxa <- pankaj2.taxa %>%
  tibble::column_to_rownames("Name")
#Remove any remaining host contamination
pankaj2.taxa <- pankaj2.taxa[which(pankaj2.taxa$Phylum != "Cyanobacteria" & pankaj2.taxa$Family != "Mitochondria"),]
pankaj2.taxa <- as.matrix(pankaj2.taxa)

#Load Oxygen Metadata
pankaj2.sort <- pankaj2[,order(colnames(pankaj2))]
pankaj2.sort <- pankaj2.sort[row.names(pankaj2.sort) %in% row.names(pankaj2.taxa),]
md.pankaj2 <- roughmetadata_pankaj2[which(unlist(roughmetadata_pankaj2$Name) %in% colnames(pankaj2.sort)),]

md.pankaj2.sort <- md.pankaj2[order(unlist(md.pankaj2$Name)),]
summary(colSums(pankaj2.sort[,-1]))

###Remove samples with read depth under 1000###
pankaj2.sort2 <- pankaj2.sort[,-c(1,which(colSums(pankaj2.sort[,-1])<1000))]
md.pankaj2.sort2 <- md.pankaj2.sort[which(md.pankaj2$Name %in% colnames(pankaj2.sort2)),]

md.pankaj2.sort2 <- md.pankaj2.sort2 %>%
  tibble::column_to_rownames("Name")

###Removes OTUs not found in any sample###
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

pankaj2.nsampls <- apply(pankaj2.sort2, 1, gt0)
pankaj2.clean <- pankaj2.sort2[which(pankaj2.nsampls>1),]

###Create phyloseq object###
OTU = otu_table(pankaj2.clean, taxa_are_rows = TRUE)
TAX = tax_table(pankaj2.taxa)
samples = sample_data(md.pankaj2.sort2)

exp5_pankaj2 <- phyloseq(OTU, TAX, samples)

###Aggregate to Genus level###
exp5_pankaj2_gen <- tax_glom(exp5_pankaj2, taxrank = 'Genus')
exp5_pankaj2_gen_prev <- filter_taxa(exp5_pankaj2_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)

########################## Figure 1H ##########################

#Get diversity indices
pankaj2.sort.flip <- t(pankaj2.sort)
simp <- vegan::diversity(pankaj2.sort.flip, "simpson")
shan <- vegan::diversity(pankaj2.sort.flip, "shannon")
indices <- as.data.frame(cbind(simp, shan))
colnames(indices) <- c("Simpson","Shannon")

#Wilcox test for HO vs NO
wilcox_alpha_pankaj2 <- t(sapply(indices, function(x) unlist(wilcox.test(x~md.pankaj2.sort$Oxygen)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_pankaj2

richness_est_pankaj2 <- indices %>%
  mutate(
    Oxygen = md.pankaj2.sort$Oxygen
  )

                                 
########################## Figure 1I ##########################

###Get differentially abundant genera with DESeq2###
diff_otu = phyloseq_to_deseq2(exp5_pankaj2_gen_prev, ~ Oxygen)
gm_mean = function(x, na.rm = TRUE){
  exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
}
diff_otu = estimateSizeFactors(diff_otu,geoMeans = apply(counts(diff_otu), 1, gm_mean))
diff_otu = DESeq(diff_otu, fitType="local")
deseq2results_otu = results(diff_otu, pAdjustMethod = "fdr")
deseq2results_otu = deseq2results_otu[order(deseq2results_otu$padj, na.last=NA), ]
sigtab_otu = deseq2results_otu[(deseq2results_otu$padj < 0.05), ]
sigtab_otu = cbind(as(sigtab_otu, "data.frame"), as(tax_table(exp5_pankaj2_gen_prev)[rownames(sigtab_otu), ], "matrix"))
print(sigtab_otu)

####Pull out specific OTUs###

exp5_pankaj2_rel_prev <- transform_sample_counts(exp5_pankaj2_gen_prev, function(x) x / sum(x) )
exp5_pankaj2_signif_feature2 <- subset_taxa(exp5_pankaj2_rel_prev, Genus == "Lactobacillus" | Genus == "Staphylococcus" | Genus == "Corynebacterium" | Genus == "Romboutsia")
exp5_pankaj2_signif_otu_feature2 <- as.data.frame(t(exp5_pankaj2_signif_feature2@otu_table))
exp5_pankaj2_signif_meta_feature2 <- as.data.frame(exp5_pankaj2_signif_feature2@sam_data)
colnames(exp5_pankaj2_signif_otu_feature2) <- c("Lactobacillus","Staphylococcus","Romboutsia","Corynebacterium")
exp5_pankaj2_diffabund_signif_feature2 <- exp5_pankaj2_signif_otu_feature2 %>%
  dplyr::mutate(SampleType = exp5_pankaj2_signif_meta_feature2$Oxygen)

#Save .csv and create boxplots in GraphPad


########################## Figure 1J ##########################

###Get relative abundance and pull out OTU table, taxonomy, and metadata from phyloseq object###
exp5_pankaj2_rel_prev <- transform_sample_counts(exp5_pankaj2_gen_prev, function(x) x / sum(x) )

exp5_pankaj2_otu_rel <- as.data.frame(t(exp5_pankaj2_rel_prev@otu_table))
exp5_pankaj2_tax_rel <- as.data.frame(exp5_pankaj2_rel_prev@tax_table)
exp5_pankaj2_meta_rel <- as.data.frame(exp5_pankaj2_rel_prev@sam_data)

###Get distance matrix and ordinate###
exp5_pankaj2_rel_bray = vegdist(exp5_pankaj2_otu_rel, method='bray')
exp5_pankaj2_rel_pcoa <- ape::pcoa(exp5_pankaj2_rel_bray)
exp5_pankaj2_rel_pcoa$values

###Set colors based on metadata###
factor_exp5 <- as.factor(exp5_pankaj2_meta_rel$Oxygen)
type_exp5 <- as.numeric(factor_exp5)
pca_colors <- c("#FB0207","#EF94A2")

###Plot###
tiff("Figure 1J.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.8, 0.4), c(-0.6, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp5_pankaj2_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp5], lwd = 1)
ordiellipse(exp5_pankaj2_rel_pcoa$vectors[,1:2], factor_exp5, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp5_pankaj2_rel_pcoa$vectors[,1:2], factor_exp5, label = TRUE)
dev.off()

###PERMANOVA and PERMDISP###
permanova_pcoa_df <- data.frame(exp5_pankaj2_meta_rel)
set.seed(1312)
permanova_pankaj2_pcoa <- vegan::adonis2(exp5_pankaj2_otu_rel ~ Oxygen, data = permanova_pcoa_df, method="bray", permutations = 10000)

pankaj2_pcoa_dist <- vegdist(exp5_pankaj2_otu_rel, method = "bray")
disp_pankaj2_pcoa <- betadisper(pankaj2_pcoa_dist, permanova_pcoa_df$Oxygen)
set.seed(1312)
permdisp_pankaj2_pcoa <- permutest(disp_pankaj2_pcoa, permutations = 10000)

print(permanova_pankaj2_pcoa)
print(permdisp_pankaj2_pcoa)


########################## Figure 1K ##########################

exp5_pankaj2_otu <- as.data.frame(t(exp5_pankaj2_gen_prev@otu_table))
exp5_pankaj2_tax <- as.data.frame(exp5_pankaj2_gen_prev@tax_table)
exp5_pankaj2_meta <- as.data.frame(exp5_pankaj2_gen_prev@sam_data)

###Hellinger transform and ordinate###
exp5_pankaj2_otu_hel <- decostand(exp5_pankaj2_otu, "hellinger")
exp5_pankaj2_otu_pca <- rda(exp5_pankaj2_otu_hel)
summary(exp5_pankaj2_otu_pca)$cont

###Set colors based on metadata###
factor_exp5 <- as.factor(exp5_pankaj2_meta$Oxygen)
type_exp5 <- as.numeric(factor_exp5)
pca_colors <- c("#FB0207","#EF94A2")

###Set priority for which labels show###
priority_exp5_pankaj2 <- colSums(exp5_pankaj2_otu)
labels_exp5_pankaj2 <- orditorp(exp5_pankaj2_otu_pca, "sp", label = exp5_pankaj2_tax$Genus, priority=priority_exp5_pankaj2)

###Plot Loading Plot###
dpi = 600
tiff("Figure 1K.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp5_pankaj2_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (44.94% Explained)", ylab="PC2 (28.10% Explained)", ylim = c(-0.6,0.3),xlim = c(-1.0,0.5))
orditorp(exp5_pankaj2_otu_pca, "sp", label = exp5_pankaj2_tax$Genus, priority=priority_exp5_pankaj2, select = (labels_exp5_pankaj2 == TRUE), cex = 0.7)
dev.off()

###Plot PCA###
tiff("Figure 1K PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.9), c(-1.0, 0.7), font = 2, font.lab = 2, xlab="PC1 (44.94% Explained)", ylab="PC2 (28.10% Explained)", type="n")
points(exp5_pankaj2_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[type_exp5], lwd = 1)
ordiellipse(exp5_pankaj2_otu_pca, factor_exp5, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp5_pankaj2_otu_pca, factor_exp5, label = TRUE)
dev.off()

###PERMANOVA and PERMDISP###
permanova_pankaj2_df <- data.frame(exp5_pankaj2_meta)
set.seed(1312)
permanova_pankaj2 <- vegan::adonis2(exp5_pankaj2_otu_hel ~ Oxygen, data = permanova_pankaj2_df, method="euclidean", permutations = 10000)

pankaj2_hel_dist <- vegdist(exp5_pankaj2_otu_hel, method = "euclidean")
disp_pankaj2 <- betadisper(pankaj2_hel_dist, permanova_pankaj2_df$Oxygen)
set.seed(1312)
permdisp_pankaj2 <- permutest(disp_pankaj2, permutations = 10000)

print(permanova_pankaj2)
print(permdisp_pankaj2)
