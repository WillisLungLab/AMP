library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(Maaslin2)
library(ecole)

pca_colors <- c("#FB0207","#B3B3B3")
deseq2_colors <- c("#FB0207","#B3B3B3")

########Import and sort OTU table########

lyso <- read_csv("exp2_lysozyme_OTU.csv")
lyso <- lyso %>%
  tibble::column_to_rownames("Name")

lyso.taxa <- read_csv("exp2_lysozyme_taxa.csv")
lyso.taxa <- lyso.taxa %>%
  tibble::column_to_rownames("Name")
lyso.taxa <- as.matrix(lyso.taxa)

roughmetadata_lyso <- read_csv("exp2_lysozyme_metadata.csv")
lyso.sort <- lyso[,order(colnames(lyso))]
md.lyso <- roughmetadata_lyso[which(unlist(roughmetadata_lyso$Name) %in% colnames(lyso.sort)),]
####Changes room air in the metadata to normoxia####
md.lyso <- data.frame(lapply(md.lyso, function(x) {gsub("RA", "NO", x)}))
md.lyso <- tibble(md.lyso)

md.lyso.sort <- md.lyso[order(unlist(md.lyso$Name)),]
####Check read depth####
summary(colSums(lyso.sort[,-1]))

########Alpha Diversity########

OTU_rough = otu_table(lyso.sort, taxa_are_rows = TRUE)

TAX_rough = tax_table(lyso.taxa)

md.lyso.sort <- md.lyso.sort %>%
  tibble::column_to_rownames("Name")
samples_rough = sample_data(md.lyso.sort)

exp5_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)
exp5_rough_gen <- tax_glom(exp5_rough, taxrank = "Genus")

fr_lyso <- prune_taxa(taxa_sums(exp5_rough_gen) > 0, exp5_rough_gen)
richness_est_lyso <- estimate_richness(fr_lyso, measures = c("Chao1", "Shannon"))
kruskal_alpha_lyso <- t(sapply(richness_est_lyso, function(x) unlist(kruskal.test(x~sample_data(fr_lyso)$SampleType)[c("estimate","p.value","statistic","conf.int")])))
kruskal_alpha_lyso

richness_est_lyso <- richness_est_lyso %>%
  mutate(
    SampleType = md.lyso.sort$SampleType
  )


plot_richness(exp5_rough_gen, x = "SampleType", measures = c("Chao1","Shannon"))
write.csv(richness_est_lyso, "2023_01_13_exp5_lyso_alpha_diversity_genus.csv")

########Processing OTU table########

####Filter out low read depth (under 1000 reads/sample)####

lyso.sort2 <- lyso.sort[,-c(1,which(colSums(lyso.sort[,-1])<1000))]
md.lyso.sort2 <- md.lyso.sort[which(md.lyso$Name %in% colnames(lyso.sort2)),]

####Filter out OTUs present in 0 samples####

gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

lyso.nsampls <- apply(lyso.sort2, 1, gt0)
lyso.clean <- lyso.sort2[which(lyso.nsampls>1),]

########Create phyloseq object########

OTU = otu_table(lyso.clean, taxa_are_rows = TRUE)
TAX = tax_table(lyso.taxa)
samples = sample_data(md.lyso.sort2)

####Aggregate to Genus level and filter out low prevalence (under 10% of samples) OTUs####
exp5_lyso_gen <- tax_glom(exp5_lyso, taxrank = 'Genus')
exp5_lyso_gen_prev <- filter_taxa(exp5_lyso_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)


######## PCA ########

###Genus###
exp5_lyso_gen_rel_prev <- transform_sample_counts(exp5_lyso_gen_prev, function(x) x / sum(x) )

exp5_lyso_otu <- as.data.frame(t(exp5_lyso_gen_rel_prev@otu_table))
exp5_lyso_tax <- as.data.frame(exp5_lyso_gen_rel_prev@tax_table)
exp5_lyso_meta <- as.data.frame(exp5_lyso_gen_rel_prev@sam_data)

exp5_lyso_otu_hel <- decostand(exp5_lyso_otu, "hellinger")
exp5_lyso_otu_pca <- rda(exp5_lyso_otu_hel)
summary(exp5_lyso_otu_pca)$cont

####Set colors for points and ellipse####
factor_exp5 <- as.factor(exp5_lyso_meta$SampleType)
factor_exp5 <- ordered(factor_exp5, levels = c("HO_LYS","HO_VEH","NO_LYS", "NO_VEH"))
factor_exp5_otherway <- ordered(factor_exp5, levels = c("NO_LYS"))
ellipse_colors <- c("#FB0207","#B3B3B3","#FB0207","#B3B3B3")
factor_exp5_ellipse <- ordered(factor_exp5, levels = c("HO_LYS","HO_VEH","NO_LYS", "NO_VEH"))

####Sets which taxa labels should be shown and which should be hidden, to prevent overcrowding####
priority_exp5_lyso <- colSums(exp5_lyso_otu_hel)
labels_exp5_lyso <- orditorp(exp5_lyso_otu_pca, "sp", label = exp5_lyso_tax$Genus, priority=priority_exp5_lyso)
dpi = 600

####PCA taxa biplot###
tiff("2022_01_13_betadiv_pca_lyso_exp5_taxa.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp5_lyso_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (36.46% Explained)", ylab="PC2 (27.87% Explained)", ylim = c(-0.5,0.8),xlim = c(-0.3,0.7))
orditorp(exp5_lyso_otu_pca, "sp", label = exp5_lyso_tax$Genus, priority=priority_exp5_lyso, select = (labels_exp5_lyso == TRUE), cex = 0.7)
dev.off()

####PCA sample plot####
tiff("2023_01_13_betadiv_pca_lyso_exp5_lyso_genus.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.8, 0.4), font = 2, font.lab = 2, xlab="PC1 (36.46% Explained)", ylab="PC2 (27.87% Explained)", type="n")
points(exp5_lyso_otu_pca, pch = 21, cex = 1.3, bg = pca_colors[factor_exp5], lwd = 1)
points(exp5_lyso_otu_pca, pch = 21, cex = 1.3, col = pca_colors[factor_exp5_otherway], lwd = 1)
ordiellipse(exp5_lyso_otu_pca, factor_exp5, col = ellipse_colors)
ordispider(exp5_lyso_otu_pca, factor_exp5, label = TRUE)
dev.off()

####PERMANOVA based on Euclidean distances####
permanova_lyso_df <- data.frame(exp5_lyso_meta)
set.seed(420)
permanova_lyso <- vegan::adonis2(exp5_lyso_otu_hel ~ SampleType, data = permanova_lyso_df, method="euclidean", permutations = 10000)
pairwise_permanova_lyso <- permanova_pairwise(exp5_lyso_otu_hel, grp = permanova_lyso_df$SampleType, permutations = 10000, method = "euclidean", padj = "fdr")

####PERMDISP based on Euclidean distances####
lyso_hel_dist <- vegdist(exp5_lyso_otu_hel, method = "euclidean")
disp_lyso <- betadisper(lyso_hel_dist, permanova_lyso_df$SampleType)
set.seed(420)
permdisp_lyso <- permutest(disp_lyso, permutations = 10000)

print(permanova_lyso)
print(pairwise_permanova_lyso)
print(permdisp_lyso)


########DESeq2########

###Genus, SampleType###

diff = phyloseq_to_deseq2(exp5_lyso_gen_prev, ~ SampleType)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")

###Creates a comparison for every combination of the 4 categories###
deseq2resultslyso1 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_VEH","HO_LYS"))
deseq2resultslyso2 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","HO_VEH","HO_LYS"))
deseq2resultslyso3 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_LYS","HO_LYS"))
deseq2resultslyso4 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_VEH","RA_LYS"))
deseq2resultslyso5 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","HO_VEH","RA_LYS"))
deseq2resultslyso6 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_VEH","HO_VEH"))

###Creates a data frame with the signficantly differentially abundant genera###
diffabund_deseq2 <- function(deseq2results){
  deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
  deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(deseq2results), ], "matrix"))
  sigtab = deseq2results[(deseq2results$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(sigtab), ], "matrix"))
  sigtabgenus = subset(sigtab, !is.na(Genus))
  print(sigtabgenus)
}
deseq2lyso1 <- diffabund_deseq2(deseq2resultslyso1)
write.csv(deseq2lyso1, "2022_01_19_exp5_lyso_diffabund_deseq_raveh_holys_stats.csv")
deseq2lyso2 <- diffabund_deseq2(deseq2resultslyso2)
write.csv(deseq2lyso2, "2022_01_19_exp5_lyso_diffabund_deseq_hoveh_holys_stats.csv")
deseq2lyso3 <- diffabund_deseq2(deseq2resultslyso3)
write.csv(deseq2lyso3, "2022_01_19_exp5_lyso_diffabund_deseq_ralys_holys_stats.csv")
deseq2lyso4 <- diffabund_deseq2(deseq2resultslyso4)
write.csv(deseq2lyso4, "2022_01_19_exp5_lyso_diffabund_deseq_raveh_ralys_stats.csv")
deseq2lyso5 <- diffabund_deseq2(deseq2resultslyso5)
write.csv(deseq2lyso5, "2022_01_19_exp5_lyso_diffabund_deseq_hoveh_ralys_stats.csv")
deseq2lyso6 <- diffabund_deseq2(deseq2resultslyso6)
write.csv(deseq2lyso6, "2022_01_19_exp5_lyso_diffabund_deseq_raveh_hoveh_stats.csv")

###Plots results###
deseq2plot <- function(sigtabgenus){
  diffabund <- ggplot(sigtabgenus,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
    geom_col() +
    coord_flip()
  diffabund2 <- diffabund + theme_bw()
  diffabund2
}
deseq2plot(deseq2lyso1)
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_raveh_holys.tiff", width = 5, height = 10, device='tiff', dpi=600)
deseq2plot(deseq2lyso2)
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_hoveh_holys.tiff", width = 5, height = 10, device='tiff', dpi=600)
deseq2plot(deseq2lyso3)
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_ralys_holys.tiff", width = 5, height = 5, device='tiff', dpi=600)
deseq2plot(deseq2lyso4)
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_raveh_ralys.tiff", width = 5, height = 10, device='tiff', dpi=600)
deseq2plot(deseq2lyso5)
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_hoveh_ralys.tiff", width = 5, height = 10, device='tiff', dpi=600)
deseq2plot(deseq2lyso6)
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_raveh_hoveh.tiff", width = 5, height = 5, device='tiff', dpi=600)


###Genus, Hyperoxia###

diff = phyloseq_to_deseq2(exp5_lyso_gen_prev, ~ Hyperoxia)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")
deseq2results = results(diff, pAdjustMethod = "fdr")
deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
sigtab = deseq2results[(deseq2results$padj < 0.05), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(sigtab), ], "matrix"))
sigtabgenus = subset(sigtab, !is.na(Genus))
print(sigtabgenus)

diffabund <- ggplot(sigtabgenus,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
  geom_col() +
  coord_flip()+
  scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
diffabund2 <- diffabund + theme_bw()
diffabund2
ggsave("2023_01_19_exp5_lyso_diffabund_deseq2_genus_hyperoxia.tiff", width = 5, height = 5, device='tiff', dpi=600)


###Make Table ###
exp5_lyso_signif_feature <- subset_taxa(exp5_lyso_gen_rel_prev, rownames(tax_table(exp5_lyso_gen_rel_prev)) %in% rownames(sigtab))
exp5_lyso_signif_otu_feature <- as.data.frame(t(exp5_lyso_signif_feature@otu_table))
exp5_lyso_signif_meta_feature <- as.data.frame(exp5_lyso_signif_feature@sam_data)
exp5_lyso_diffabund_signif_feature <- exp5_lyso_signif_otu_feature %>%
  dplyr::mutate(SampleType = exp5_lyso_signif_meta_feature$SampleType)

write.csv(exp5_lyso_diffabund_signif_feature, "2022_01_12_exp5_lyso_diffabund_deseq_genus_hyperoxia.csv")
write.csv(sigtab, "2022_01_12_exp5_lyso_diffabund_deseq_genus_hyperoxia_stats.csv")

########MaAslin2########

####Process data to run MaAsLin2 at genus level####

lyso.clean.taxa <- as.data.frame(lyso.taxa[which(unlist(rownames(lyso.taxa)) %in% rownames(lyso.clean)),])
lyso.clean.taxa[,6][is.na(lyso.clean.taxa[,6])] = "Unidentified"

lyso.clean.genus <- lyso.clean %>%
  dplyr::mutate(Genus = lyso.clean.taxa$Genus) %>%
  group_by(Genus) %>% 
  summarize(across(.cols = everything(), .fns = sum)) %>%
  tibble::column_to_rownames("Genus")

###Genus, SampleType###

lyso.clean.genus.flip <- t(lyso.clean.genus)

fit_data_exp5_lyso = Maaslin2(
  input_data = lyso.clean.genus.flip,
  input_metadata = md.lyso.sort2,
  output = "exp5_lyso_output_genus",
  max_significance = 0.05,
  fixed_effects = c("Hyperoxia","Lyso")
)

###Abundance Boxplot###


exp5_lyso_rel_prev <- transform_sample_counts(exp5_lyso_prev, function(x) x / sum(x) )
exp5_lyso_gen_rel_prev <- transform_sample_counts(exp5_lyso_gen_prev, function(x) x / sum(x) )

exp5_lyso_diffabund_subset <- subset_taxa(exp5_lyso_gen_rel_prev, rownames(tax_table(exp5_lyso_gen_rel_prev)) %in% rownames(sigtabgenus))

exp5_lyso_diffabund_subset <- subset_taxa(exp5_lyso_gen_rel_prev, rownames(tax_table(exp5_lyso_gen_rel_prev)) %in% c("OTU679"))
exp5_lyso_diffabund_df <- psmelt(exp5_lyso_diffabund_subset)

diffabund_boxplot <- ggplot(data = exp5_lyso_diffabund_df, aes(x = SampleType, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = SampleType), height = 0, width = .2) +
  labs(x = "", y = "Abundance" ) +
  facet_wrap(~ Genus, scales = "free") +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))

diffabund_boxplot + scale_color_manual(values = c("#FB0207","#FB0207","#B3B3B3","#B3B3B3"))
ggsave("2023_01_12_exp5_1_diffabund_genus_boxplot.tiff", width = 10, height = 5, device='tiff', dpi=600)

exp5_lyso_signif <- subset_taxa(exp5_lyso_gen_rel_prev, rownames(tax_table(exp5_lyso_gen_rel_prev)) %in% c("OTU4893"))
exp5_lyso_signif_otu <- as.data.frame(t(exp5_lyso_signif@otu_table))
exp5_lyso_signif_meta <- as.data.frame(exp5_lyso_signif@sam_data)
exp5_lyso_signif_otu <- exp5_lyso_signif_otu %>%
  dplyr::mutate(SampleType = exp5_lyso_signif_meta$SampleType)
write.csv(exp5_lyso_signif_otu, "2023_01_20_exp5_lyso_genus_streptococcus.csv")

deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(deseq2results), ], "matrix"))
deseq2results_df <- deseq2results_df[c(1,9,30,39,4,136),]
write.csv(deseq2results_df, "2023_01_13_exp5_lyso_genus_signif_stats.csv")



