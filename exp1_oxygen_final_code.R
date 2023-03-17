library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(mvabund)
library(Maaslin2)

###Set Colors###
betadiv_colors <- c("#EF94A2", "#FB0207")
deseq2_colors <- c("#FB0207","#EF94A2")

###Load Data and Arrange###
pankaj2 <- read_csv("exp1_oxygen_OTU.csv")
pankaj2 <- pankaj2 %>%
  tibble::column_to_rownames("Name")

pankaj2.taxa <- read_csv("exp1_oxygen_taxa.csv")
pankaj2.taxa <- pankaj2.taxa %>%
  tibble::column_to_rownames("Name")
pankaj2.taxa <- as.matrix(pankaj2.taxa)

roughmetadata_pankaj2 <- read_table("exp1_oxygen_metadata.txt")

pankaj2.sort <- pankaj2[,order(colnames(pankaj2))]
md.pankaj2 <- roughmetadata_pankaj2[which(unlist(roughmetadata_pankaj2$Name) %in% colnames(pankaj2.sort)),]

md.pankaj2.sort <- md.pankaj2[order(unlist(md.pankaj2$Name)),]
summary(colSums(pankaj2.sort[,-1]))


###Alpha Diversity###

OTU_rough = otu_table(pankaj2.sort, taxa_are_rows = TRUE)

TAX_rough = tax_table(pankaj2.taxa)

md.pankaj2.sort <- md.pankaj2.sort %>%
  tibble::column_to_rownames("Name")
md.pankaj2.sort$Oxygen <- as.factor(md.pankaj2.sort$Oxygen)
samples_rough = sample_data(md.pankaj2.sort)

exp5_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)
exp5_rough_gen <- tax_glom(exp5_rough, taxrank = "Genus")

###Genus###

fr_pankaj2 <- prune_taxa(taxa_sums(exp5_rough_gen) > 0, exp5_rough_gen)
richness_est_pankaj2 <- estimate_richness(fr_pankaj2, measures = c("Chao1", "Shannon"))
wilcox_alpha_pankaj2 <- t(sapply(richness_est_pankaj2, function(x) unlist(wilcox.test(x~sample_data(fr_pankaj2)$Oxygen)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_pankaj2
richness_est_pankaj2 <- richness_est_pankaj2 %>%
  mutate(
    Oxygen = md.pankaj2.sort$Oxygen
  )
plot_richness(exp5_rough_gen, x = "Oxygen", measures = c("Chao1","Shannon")) + theme_bw()
write.csv(richness_est_pankaj2, "Pankaj2 Alpha Diversity.csv")

###Remove samples with read depth under 1000###
pankaj2.sort2 <- pankaj2.sort[,-c(1,which(colSums(pankaj2.sort[,-1])<1000))]
md.pankaj2.sort2 <- md.pankaj2.sort[which(md.pankaj2$Name %in% colnames(pankaj2.sort2)),]

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

###Filter to Genus level and filter out OTUs present in <10% of samples###
exp5_pankaj2_gen <- tax_glom(exp5_pankaj2, taxrank = 'Genus')
exp5_pankaj2_gen_prev <- filter_taxa(exp5_pankaj2_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)


###Prepare PCA###

exp5_pankaj2_otu <- as.data.frame(t(exp5_pankaj2_gen_prev@otu_table))
exp5_pankaj2_tax <- as.data.frame(exp5_pankaj2_gen_prev@tax_table)
exp5_pankaj2_meta <- as.data.frame(exp5_pankaj2_gen_prev@sam_data)

exp5_pankaj2_otu_hel <- decostand(exp5_pankaj2_otu, "hellinger")
exp5_pankaj2_otu_pca <- rda(exp5_pankaj2_otu_hel)
summary(exp5_pankaj2_otu_pca)$cont

factor_exp5 <- as.factor(exp5_pankaj2_meta$Oxygen)
type_exp5 <- as.numeric(factor_exp5)
pca_colors <- c("#FB0207","#EF94A2")

priority_exp5_pankaj2 <- colSums(exp5_pankaj2_otu)
labels_exp5_pankaj2 <- orditorp(exp5_pankaj2_otu_pca, "sp", label = exp5_pankaj2_tax$Genus, priority=priority_exp5_pankaj2)
dpi = 600
tiff("PCA Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp5_pankaj2_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (44.64% Explained)", ylab="PC2 (27.97% Explained)", ylim = c(-0.8,0.3),xlim = c(-1.0,0.5))
orditorp(exp5_pankaj2_otu_pca, "sp", label = exp5_pankaj2_tax$Genus, priority=priority_exp5_pankaj2, select = (labels_exp5_pankaj2 == TRUE), cex = 0.7)
dev.off()

tiff("Beta Diversity PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.4, 0.9), c(-0.9, 0.7), font = 2, font.lab = 2, xlab="PC1 (44.64% Explained)", ylab="PC2 (27.97% Explained)", type="n")
points(exp5_pankaj2_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[type_exp5], lwd = 1)
ordispider(exp5_pankaj2_otu_pca, factor_exp5, label = TRUE)
ordiellipse(exp5_pankaj2_otu_pca, factor_exp5, col = pca_colors)
dev.off()

permanova_pankaj2_df <- data.frame(exp5_pankaj2_meta)
set.seed(1312)
permanova_pankaj2 <- vegan::adonis2(exp5_pankaj2_otu_hel ~ Oxygen, data = permanova_pankaj2_df, method="euclidean", permutations = 10000)

pankaj2_hel_dist <- vegdist(exp5_pankaj2_otu_hel, method = "euclidean")
disp_pankaj2 <- betadisper(pankaj2_hel_dist, permanova_pankaj2_df$Oxygen)
set.seed(1312)
permdisp_pankaj2 <- permutest(disp_pankaj2, permutations = 10000)

print(permanova_pankaj2)
print(permdisp_pankaj2)

###PCoA###
exp5_pankaj2_rel_prev <- transform_sample_counts(exp5_pankaj2_gen_prev, function(x) x / sum(x) )

exp5_pankaj2_otu_rel <- as.data.frame(t(exp5_pankaj2_rel_prev@otu_table))
exp5_pankaj2_tax_rel <- as.data.frame(exp5_pankaj2_rel_prev@tax_table)
exp5_pankaj2_meta_rel <- as.data.frame(exp5_pankaj2_rel_prev@sam_data)

exp5_pankaj2_rel_bray = vegdist(exp5_pankaj2_otu_rel, method='bray')
exp5_pankaj2_rel_pcoa <- ape::pcoa(exp5_pankaj2_rel_bray)
exp5_pankaj2_rel_pcoa$values

factor_exp5 <- as.factor(exp5_pankaj2_meta$Oxygen)
type_exp5 <- as.numeric(factor_exp5)
pca_colors <- c("#FB0207","#EF94A2")

tiff("Beta Diversity PCoA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.6, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp5_pankaj2_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp5], lwd = 1)
ordispider(exp5_pankaj2_rel_pcoa$vectors[,1:2], factor_exp5, label = TRUE)
ordiellipse(exp5_pankaj2_rel_pcoa$vectors[,1:2], factor_exp5, col = pca_colors)
dev.off()

permanova_pcoa_df <- data.frame(exp5_pankaj2_meta_rel)
set.seed(1312)
permanova_pankaj2_pcoa <- vegan::adonis2(exp5_pankaj2_otu_rel ~ Oxygen, data = permanova_pcoa_df, method="bray", permutations = 10000)

pankaj2_pcoa_dist <- vegdist(exp5_pankaj2_otu_rel, method = "bray")
disp_pankaj2_pcoa <- betadisper(pankaj2_hel_dist, permanova_pcoa_df$Oxygen)
set.seed(1312)
permdisp_pankaj2_pcoa <- permutest(disp_pankaj2_pcoa, permutations = 10000)

print(permanova_pankaj2)
print(permdisp_pankaj2)

###DESeq2###

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
write.csv(sigtab_otu,"DESeq2 Signficant.csv")

###Table###

pankaj2.clean.taxa <- as.data.frame(pankaj2.taxa[which(unlist(rownames(pankaj2.taxa)) %in% rownames(pankaj2.clean)),])
pankaj2.clean.taxa[,6][is.na(pankaj2.clean.taxa[,6])] = "Unidentified"

pankaj2.clean.genus <- pankaj2.clean %>%
  dplyr::mutate(Genus = pankaj2.clean.taxa$Genus) %>%
  group_by(Genus) %>% 
  summarize(across(.cols = everything(), .fns = sum)) %>%
  tibble::column_to_rownames("Genus")

###MaAslin2###
pankaj2.clean.taxa <- as.data.frame(pankaj2.taxa[which(unlist(rownames(pankaj2.taxa)) %in% rownames(pankaj2.clean)),])
pankaj2.clean.taxa[,6][is.na(pankaj2.clean.taxa[,6])] = "Unidentified"

pankaj2.clean.genus <- pankaj2.clean %>%
  dplyr::mutate(Genus = pankaj2.clean.taxa$Genus) %>%
  group_by(Genus) %>% 
  summarize(across(.cols = everything(), .fns = sum)) %>%
  tibble::column_to_rownames("Genus")

pankaj2.clean.flip <- t(pankaj2.clean.genus)

fit_data_exp5_pankaj2 = Maaslin2(
  input_data = pankaj2.clean.flip,
  input_metadata = md.pankaj2.sort2,
  output = "exp5_pankaj2_output",
  max_significance = 0.05,
  fixed_effects = "Oxygen"
)

###Write .csv file to create boxplot in GraphPad###
exp5_pankaj2_rel_prev <- transform_sample_counts(exp5_pankaj2_gen_prev, function(x) x / sum(x) )

exp5_pankaj2_signif_feature <- subset_taxa(exp5_pankaj2_rel_prev, rownames(tax_table(exp5_pankaj2_rel_prev)) %in% rownames(sigtab_otu))
exp5_pankaj2_signif_otu_feature <- as.data.frame(t(exp5_pankaj2_signif_feature@otu_table))
exp5_pankaj2_signif_meta_feature <- as.data.frame(exp5_pankaj2_signif_feature@sam_data)
colnames(exp5_pankaj2_signif_otu_feature) <- c("OTU002: Staphylococcus","OTU013: Romboutsia","OTU024: Proteus","OTU026: Bacteroides","OTU077: Herbaspirillum","OTU124: Corynebacterium")
exp5_pankaj2_diffabund_signif_feature <- exp5_pankaj2_signif_otu_feature %>%
  dplyr::mutate(SampleType = exp5_pankaj2_signif_meta_feature$Oxygen)

exp5_pankaj_diffabund_subset <- subset_taxa(exp5_pankaj2_rel_prev, rownames(tax_table(exp5_pankaj2_rel_prev)) %in% rownames(sigtab_otu))
exp5_pankaj_diffabund_df <- psmelt(exp5_pankaj_diffabund_subset)

write.csv(exp5_pankaj2_diffabund_signif_feature, "DESeq2 Significant.csv")

###Preliminary Boxplot###

diffabund_boxplot <- ggplot(data = exp5_pankaj_diffabund_df, aes(x = Oxygen, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = Oxygen), height = 0, width = .2) +
  labs(x = "", y = "Abundance" ) +
  facet_wrap(~ Genus, scales = "free") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

diffabund_boxplot + scale_color_manual(values = c("#FB0207","#EF94A2"))
