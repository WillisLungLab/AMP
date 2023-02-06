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
pankaj <- read_csv("exp1_oxygen_OTU.csv")
pankaj <- pankaj %>%
  tibble::column_to_rownames("Name")

pankaj.taxa <- read_csv("exp1_oxygen_taxa.csv")
pankaj.taxa <- pankaj.taxa %>%
  tibble::column_to_rownames("Name")
pankaj.taxa <- as.matrix(pankaj.taxa)

roughmetadata_pankaj <- read_table("exp1_oxygen_metadata.txt")

pankaj.sort <- pankaj[,order(colnames(pankaj))]
md.pankaj <- roughmetadata_pankaj[which(unlist(roughmetadata_pankaj$Name) %in% colnames(pankaj.sort)),]

md.pankaj.sort <- md.pankaj[order(unlist(md.pankaj$Name)),]
summary(colSums(pankaj.sort[,-1]))

### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. ###
### 6432   13803   15791   16233   20112   22219 ###

###Alpha Diversity###

OTU_rough = otu_table(pankaj.sort, taxa_are_rows = TRUE)

TAX_rough = tax_table(pankaj.taxa)

md.pankaj.sort <- md.pankaj.sort %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.pankaj.sort)

exp5_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)
exp5_rough_gen <- tax_glom(exp5_rough, taxrank = "Genus")

###Genus###
fr_pankaj <- prune_taxa(taxa_sums(exp5_rough_gen) > 0, exp5_rough_gen)
richness_est_pankaj <- estimate_richness(fr_pankaj, measures = c("Chao1", "Shannon"))
wilcox_alpha_pankaj <- t(sapply(richness_est_pankaj, function(x) unlist(wilcox.test(x~sample_data(fr_pankaj)$SampleType)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_pankaj
richness_est_pankaj <- richness_est_pankaj %>%
  mutate(
    SampleType = md.pankaj.sort$SampleType
  )
plot_richness(exp5_rough_gen, x = "SampleType", measures = c("Chao1","Shannon"))
write.csv(richness_est_pankaj, "2023_01_12_exp5_pankaj_alpha_diversity_genus.csv")

###Remove samples with read depth under 1000###
pankaj.sort2 <- pankaj.sort[,-c(1,which(colSums(pankaj.sort[,-1])<1000))]
md.pankaj.sort2 <- md.pankaj.sort[which(md.pankaj$Name %in% colnames(pankaj.sort2)),]

###Removes OTUs not found in any sample###
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

pankaj.nsampls <- apply(pankaj.sort2, 1, gt0)
pankaj.clean <- pankaj.sort2[which(pankaj.nsampls>1),]

###Create phyloseq object###
OTU = otu_table(pankaj.clean, taxa_are_rows = TRUE)
TAX = tax_table(pankaj.taxa)
samples = sample_data(md.pankaj.sort2)

exp5_pankaj <- phyloseq(OTU, TAX, samples)

###Filter to Genus level and filter out OTUs present in <10% of samples###
exp5_pankaj_gen <- tax_glom(exp5_pankaj, taxrank = 'Genus')
exp5_pankaj_gen_prev <- filter_taxa(exp5_pankaj_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)


###Prepare PCA###

exp5_pankaj_rel_prev <- transform_sample_counts(exp5_pankaj_prev, function(x) x / sum(x) )
exp5_pankaj_gen_rel_prev <- transform_sample_counts(exp5_pankaj_gen_prev, function(x) x / sum(x) )

exp5_pankaj_otu <- as.data.frame(t(exp5_pankaj_gen_rel_prev@otu_table))
exp5_pankaj_tax <- as.data.frame(exp5_pankaj_gen_rel_prev@tax_table)
exp5_pankaj_meta <- as.data.frame(exp5_pankaj_gen_rel_prev@sam_data)

exp5_pankaj_otu_hel <- decostand(exp5_pankaj_otu, "hellinger")
exp5_pankaj_otu_pca <- rda(exp5_pankaj_otu_hel)
summary(exp5_pankaj_otu_pca)$cont

factor_exp5 <- as.factor(exp5_pankaj_meta$SampleType)
type_exp5 <- as.numeric(factor_exp5)
pca_colors <- c("#FB0207","#EF94A2")

priority_exp5_pankaj <- colSums(exp5_pankaj_otu_hel)
labels_exp5_pankaj <- orditorp(exp5_pankaj_otu_pca, "sp", label = exp5_pankaj_tax$Genus, priority=priority_exp5_pankaj)
dpi = 600
tiff("2022_01_12_betadiv_pca_pankaj_exp5_taxa_genus.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp5_pankaj_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (54.57% Explained)", ylab="PC2 (29.94% Explained)", ylim = c(-0.3,0.3),xlim = c(-0.5,1.2))
orditorp(exp5_pankaj_otu_pca, "sp", label = exp5_pankaj_tax$Genus, priority=priority_exp5_pankaj, select = (labels_exp5_pankaj == TRUE), cex = 0.7)
dev.off()

tiff("2022_01_12_betadiv_pca_pankaj_exp5_ho_ra_genus.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.5), c(-0.6, 0.5), font = 2, font.lab = 2, xlab="PC1 (54.57% Explained)", ylab="PC2 (29.94% Explained)", type="n")
points(exp5_pankaj_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[type_exp5], lwd = 1)
ordispider(exp5_pankaj_otu_pca, factor_exp5, label = FALSE)
ordiellipse(exp5_pankaj_otu_pca, factor_exp5, col = pca_colors)
dev.off()

permanova_pankaj_df <- data.frame(exp5_pankaj_meta)
set.seed(420)
permanova_pankaj <- vegan::adonis2(exp5_pankaj_otu_hel ~ SampleType, data = permanova_pankaj_df, method="euclidean", permutations = 10000)

pankaj_hel_dist <- vegdist(exp5_pankaj_otu_hel, method = "euclidean")
disp_pankaj <- betadisper(pankaj_hel_dist, permanova_pankaj_df$SampleType)
set.seed(420)
permdisp_pankaj <- permutest(disp_pankaj, permutations = 10000)

print(permanova_pankaj)
print(permdisp_pankaj)

###DESeq2###

diff_otu = phyloseq_to_deseq2(exp5_pankaj_prev, ~ SampleType)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans_otu = apply(counts(diff_otu), 1, gm_mean)
diff_otu = estimateSizeFactors(diff_otu, geoMeans_otu = geoMeans_otu)
diff_otu = DESeq(diff_otu, fitType="local")
deseq2results_otu = results(diff_otu, pAdjustMethod = "fdr")
deseq2results_otu = deseq2results_otu[order(deseq2results_otu$padj, na.last=NA), ]
sigtab_otu = deseq2results_otu[(deseq2results_otu$padj < 0.05), ]
sigtab_otu = cbind(as(sigtab_otu, "data.frame"), as(tax_table(exp5_pankaj)[rownames(sigtab_otu), ], "matrix"))
print(sigtabgenus_otu)
write.csv(sigtab_otu,"2022_01_12_exp5_pankaj_deseq2_otu.csv")

diffabund2
ggsave("2023_01_19_exp5_1_diffabund_deseq2.tiff", width = 5, height = 5, device='tiff', dpi=600)

###Genus Level###

diff_genus = phyloseq_to_deseq2(exp5_pankaj_gen_prev, ~ SampleType)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans_genus = apply(counts(diff), 1, gm_mean)
diff_genus = estimateSizeFactors(diff_genus, geoMeans = geoMeans)
diff_genus = DESeq(diff_genus, fitType="local")
deseq2results = results(diff_genus, pAdjustMethod = "fdr")
deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
sigtab_genus = deseq2results[(deseq2results$padj < 0.05), ]
sigtab_genus = cbind(as(sigtab_genus, "data.frame"), as(tax_table(exp5_pankaj_gen)[rownames(sigtab_genus), ], "matrix"))
print(sigtab_genus)
write.csv(sigtab_genus,"2022_01_12_exp5_pankaj_deseq2_genus.csv")

###Table###

pankaj.clean.taxa <- as.data.frame(pankaj.taxa[which(unlist(rownames(pankaj.taxa)) %in% rownames(pankaj.clean)),])
pankaj.clean.taxa[,6][is.na(pankaj.clean.taxa[,6])] = "Unidentified"

pankaj.clean.genus <- pankaj.clean %>%
  dplyr::mutate(Genus = pankaj.clean.taxa$Genus) %>%
  group_by(Genus) %>% 
  summarize(across(.cols = everything(), .fns = sum)) %>%
  tibble::column_to_rownames("Genus")

###MaAslin2###
pankaj.clean.flip <- t(pankaj.clean)

fit_data_exp5_pankaj = Maaslin2(
  input_data = pankaj.clean.flip,
  input_metadata = md.pankaj.sort2,
  output = "exp5_pankaj_output",
  max_significance = 0.05,
  fixed_effects = "SampleType"
)

###Write .csv file to create boxplot in GraphPad###
exp5_pankaj_rel_prev <- transform_sample_counts(exp5_pankaj_prev, function(x) x / sum(x) )

exp5_pankaj_signif_feature <- subset_taxa(exp5_pankaj_rel_prev, rownames(tax_table(exp5_pankaj_rel_prev)) %in% rownames(sigtab))
exp5_pankaj_signif_otu_feature <- as.data.frame(t(exp5_pankaj_signif_feature@otu_table))
exp5_pankaj_signif_meta_feature <- as.data.frame(exp5_pankaj_signif_feature@sam_data)
colnames(exp5_pankaj_signif_otu_feature) <- c("OTU013: Staphylococcus","OTU083: Staphylococcus","OTU085: Turicibacter")
exp5_pankaj_diffabund_signif_feature <- exp5_pankaj_signif_otu_feature %>%
  dplyr::mutate(SampleType = exp5_pankaj_signif_meta_feature$SampleType)

write.csv(exp5_pankaj_diffabund_signif_feature, "exp5_pankaj_diffabund_signif_feature.csv")
