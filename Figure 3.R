library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(Maaslin2)
library(ecole)

#If you used rbiom before this, detach the package now or it'll conflict with other packages
detach("package:rbiom", unload = TRUE)

pca_colors <- c("#FB0207","#B3B3B3")
deseq2_colors <- c("#FB0207","#B3B3B3")

########Import and sort OTU table########

lyso <- read_csv("Exp5_2_Decontam.csv")
lyso <- lyso %>%
  tibble::column_to_rownames("Name")

lyso.taxa <- read_csv("Exp5_2_Decontam_Taxa.csv")
lyso.taxa <- lyso.taxa %>%
  tibble::column_to_rownames("Name")
lyso.taxa <- lyso.taxa[which(lyso.taxa$Phylum != "Cyanobacteria" & lyso.taxa$Family != "Mitochondria"),]
lyso.taxa <- as.matrix(lyso.taxa)

roughmetadata_lyso <- read_csv("Exp5_2_Metadata.csv")
lyso.sort <- lyso[,order(colnames(lyso))]
lyso.sort <- lyso.sort[row.names(lyso.sort) %in% row.names(lyso.taxa),]
md.lyso <- roughmetadata_lyso[which(unlist(roughmetadata_lyso$Name) %in% colnames(lyso.sort)),]
####Changes room air in the metadata to normoxia####
md.lyso <- data.frame(lapply(md.lyso, function(x) {gsub("NO", "NO", x)}))
md.lyso <- tibble(md.lyso)

md.lyso.sort <- md.lyso[order(unlist(md.lyso$Name)),]
####Check read depth####
summary(colSums(lyso.sort[,-1]))

lyso.sort2 <- lyso.sort[,-c(1,which(colSums(lyso.sort[,-1])<1000))]

md.lyso.sort2 <- md.lyso.sort[which(md.lyso$Name %in% colnames(lyso.sort2)),]
md.lyso.sort2 <- md.lyso.sort2 %>%
  tibble::column_to_rownames("Name")

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

exp5_lyso <- phyloseq(OTU, TAX, samples)

####Aggregate to Genus level####
exp5_lyso_gen <- tax_glom(exp5_lyso, taxrank = 'Genus')
exp5_lyso_gen_prev <- filter_taxa(exp5_lyso_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)



########################## Figure 3H ###########################

lyso.sort.flip <- t(lyso.sort)
simp <- vegan::diversity(lyso.sort.flip, "simpson")
shan <- vegan::diversity(lyso.sort.flip, "shannon")
indices <- as.data.frame(cbind(simp, shan))
colnames(indices) <- c("Simpson","Shannon")

wilcox_alpha_lyso <- t(sapply(indices, function(x) unlist(kruskal.test(x~md.lyso.sort$SampleType)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_lyso

richness_est_lyso <- indices %>%
  mutate(
    SampleType = md.lyso.sort$SampleType
  )

#Save .csv and plot in GraphPad


########################## Figure 3I ###########################

exp5_lyso_signif_feature <- subset_taxa(exp5_lyso_gen_rel_prev, Genus == "Lactobacillus" | Genus == "Staphylococcus"| Genus == "Corynebacterium" | Genus == "Enterobacter" | Genus == "Erwinia" | Genus == "Klebsiella" | Genus == "Romboutsia")
exp5_lyso_signif_otu_feature <- as.data.frame(t(exp5_lyso_signif_feature@otu_table))
exp5_lyso_signif_meta_feature <- as.data.frame(exp5_lyso_signif_feature@sam_data)
colnames(exp5_lyso_signif_otu_feature) <- c("OTU001:Lactobacillus","OTU003:Staphylococcus","OTU2274:Klebsiella","OTU428:Corynebacterium","OTU573:Enterobacter","OTU966:Erwinia")
exp5_lyso_diffabund_signif_feature <- exp5_lyso_signif_otu_feature %>%
  dplyr::mutate(SampleType = exp5_lyso_signif_meta_feature$SampleType)

#Save .csv, plot in GraphPad


########################## Figure 3J ###########################

exp5_lyso_rel_prev <- transform_sample_counts(exp5_lyso_gen_prev, function(x) x / sum(x) )

exp5_lyso_otu_rel <- as.data.frame(t(exp5_lyso_rel_prev@otu_table))
exp5_lyso_tax_rel <- as.data.frame(exp5_lyso_rel_prev@tax_table)
exp5_lyso_meta_rel <- as.data.frame(exp5_lyso_rel_prev@sam_data)

###Get distance matrix and ordinate###
exp5_lyso_rel_bray = vegdist(exp5_lyso_otu_rel, method='bray')
exp5_lyso_rel_pcoa <- ape::pcoa(exp5_lyso_rel_bray)
exp5_lyso_rel_pcoa$values

###Set colors###
factor_exp5 <- as.factor(exp5_lyso_meta_rel$SampleType)
factor_exp5 <- ordered(factor_exp5, levels = c("HO_LYS","HO_VEH","NO_LYS", "NO_VEH"))
factor_exp5_otherway <- ordered(factor_exp5, levels = c("NO_LYS"))
ellipse_colors <- c("#FB0207","#B3B3B3","#FB0207","#B3B3B3")
factor_exp5_ellipse <- ordered(factor_exp5, levels = c("HO_LYS","HO_VEH","NO_LYS", "NO_VEH"))

tiff("Figure 3J.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.7, 0.3), c(-0.5, 0.5), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(exp5_lyso_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[factor_exp5], lwd = 1)
points(exp5_lyso_rel_pcoa$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors[factor_exp5_otherway], lwd = 1)
ordiellipse(exp5_lyso_rel_pcoa$vectors[,1:2], factor_exp5, kind = "se", conf = 0.95,col = pca_colors)
ordispider(exp5_lyso_rel_pcoa$vectors[,1:2], factor_exp5, label = TRUE)
dev.off()

###PERMANOVA and PERMDISP###
permanova_pcoa_df <- data.frame(exp5_lyso_meta_rel)
set.seed(1312)
permanova_lyso_pcoa <- vegan::adonis2(exp5_lyso_otu_rel ~ SampleType, data = permanova_pcoa_df, method="bray", permutations = 10000)

pairwise_permanova_lyso <- permanova_pairwise(exp5_lyso_otu_rel, grp = permanova_pcoa_df$SampleType, permutations = 10000, method = "bray", padj = "fdr")

lyso_pcoa_dist <- vegdist(exp5_lyso_otu_rel, method = "bray")
disp_lyso_pcoa <- betadisper(lyso_pcoa_dist, permanova_pcoa_df$SampleType)
set.seed(1312)
permdisp_lyso_pcoa <- permutest(disp_lyso_pcoa, permutations = 10000)

print(permanova_lyso)
print(pairwise_permanova_lyso)
print(permdisp_lyso)


########################## Figure 3K ###########################

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
tiff("Figure 3K Loading Plot.tif", width=5*dpi, height=5*dpi, res=dpi)
biplot(exp5_lyso_otu_pca, display = c("species"), type = c("points"), col = "grey",font = 2, font.lab = 2, xlab="PC1 (36.60% Explained)", ylab="PC2 (27.73% Explained)", ylim = c(-0.5,0.8),xlim = c(-0.3,0.7))
orditorp(exp5_lyso_otu_pca, "sp", label = exp5_lyso_tax$Genus, priority=priority_exp5_lyso, select = (labels_exp5_lyso == TRUE), cex = 0.7)
dev.off()

####PCA sample plot####
tiff("Figure 3K PCA.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.6, 0.4), c(-0.8, 0.4), font = 2, font.lab = 2, xlab="PC1 (36.60% Explained)", ylab="PC2 (27.73% Explained)", type="n")
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
