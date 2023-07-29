library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(mvabund)
library(Maaslin2)
library(rbiom)

########################## Figure S7A ##########################

###Renyi diversity profile using Hill numbers 1-4###
set.seed(420)
renyihill <- renyi(lyso.sort.flip, scales = c(1, 2, 3, 4),hill = TRUE)

###Add metadata factors###
renyihill2 <- renyihill %>%
  mutate(
    SampleType = md.lyso.sort$SampleType,
    Lyso = md.lyso.sort$Lyso,
    Hyperoxia = md.lyso.sort$Hyperoxia,
  )

renyihill2 <- rownames_to_column(renyihill2, var = "Name")
renyihill2_long <- renyihill2 %>% pivot_longer(cols=c(2:5),
                                               names_to='Hill',
                                               values_to='Value')

###Plot###
ggplot(renyihill2_long,aes(x=Hill,y=Value,group=Name))+
  geom_line(aes(color=Hyperoxia)) +
  geom_point(aes(color=Hyperoxia, shape=Lyso), size = 3) +
  scale_color_manual(values = c("#FB0207","#B3B3B3")) +
  scale_shape_manual(values = c(16,1)) +
  theme_bw()
#Save


########################## Figure S7B ##########################

###Set Colors###
pca_colors <- c("#FB0207","#B3B3B3")
deseq2_colors <- c("#FB0207","#EF94A2")

#Load tree compatible OTU table
lyso <- lyso %>%
  tibble::column_to_rownames("Name")

#Load associated taxonomy table
lyso.taxa <- lyso.taxa %>%
  tibble::column_to_rownames("Name")
lyso.taxa <- as.matrix(lyso.taxa)

#Load Metadata

lyso.sort <- lyso[,order(colnames(lyso))]
md.lyso <- roughmetadata_lyso[which(unlist(roughmetadata_lyso$Name) %in% colnames(lyso.sort)),]

md.lyso.sort <- md.lyso[order(unlist(md.lyso$Name)),]
summary(colSums(lyso.sort[,-1]))

md.lyso.sort <- md.lyso.sort %>%
  tibble::column_to_rownames("Name")

###Remove samples with read depth under 1000###
lyso.sort2 <- lyso.sort[,-c(1,which(colSums(lyso.sort[,-1])<1000))]
md.lyso.sort2 <- md.lyso.sort[which(md.lyso$Name %in% colnames(lyso.sort2)),]

###Removes OTUs not found in any sample###
gt0 <- function(vec){
  v <- as.numeric(vec)
  s <- sum(v>0)
  return(s)
}

lyso.nsampls <- apply(lyso.sort2, 1, gt0)
lyso.clean <- lyso.sort2[which(lyso.nsampls>1),]

#Load Tree
tree <- ape::read.nexus("AMP_All_Tree.nxs")

labels <- as.data.frame(tree[["tip.label"]])
labels[c('Taxa', 'Name')] <- str_split_fixed(labels$`tree[["tip.label"]]`, ";", 2)
newlabel <- gsub("(^')|('$)", "", labels$Name)
tree[["tip.label"]] <- newlabel

lyso.clean.new <- lyso.clean[row.names(lyso.clean) %in% newlabel,]
lyso.taxa.new <- lyso.taxa[row.names(lyso.taxa) %in% newlabel,]

###Create phyloseq object###
OTU_lyso_withtree = otu_table(lyso.clean.new, taxa_are_rows = TRUE)
TAX_lyso_withtree = tax_table(lyso.taxa.new)
samples_lyso_withtree = sample_data(md.lyso.sort2)
phylo_lyso_withtree <- phy_tree(tree)

exp5_lyso_withtree <- phyloseq(OTU_lyso_withtree, TAX_lyso_withtree, samples_lyso_withtree, phylo_lyso_withtree)
exp5_lyso_gen_withtree <- tax_glom(exp5_lyso_withtree, taxrank = 'Genus')

plot_tree(exp5_lyso_gen_withtree, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color = "Phylum")
#Save

exp5_lyso_rel_withtree <- transform_sample_counts(exp5_lyso_gen_withtree, function(x) x / sum(x) )

###Get UniFrac distance and ordinate###
uni_lyso <- UniFrac(exp5_lyso_rel_withtree, weighted = TRUE)

ord_lyso = ordinate(exp5_lyso_rel, "PCoA", "unifrac", weighted=TRUE)
betadiv_lyso <- plot_ordination(exp5_lyso_gen_lyso, ord, color="SampleType", shape="SampleType") + theme_bw()
betadiv_lyso + geom_point(size = 3) + stat_ellipse()


exp5_lyso_meta <- as.data.frame(exp5_lyso_rel@sam_data)
exp5_lyso_phylo <- phy_tree(exp5_lyso_rel)

###Set colors###
factor_exp5_lyso_withtree <- as.factor(exp5_lyso_meta$SampleType)
factor_exp5_lys <- ordered(factor_exp5_lyso_withtree, levels = c("HO_LYS","NO_LYS"))
factor_exp5_veh <- ordered(factor_exp5_lyso_withtree, levels = c("HO_VEH","NO_VEH"))
ellipse_colors_lyso_withtree <- c("#FB0207","#B3B3B3","#FB0207","#B3B3B3")
factor_exp5_ellipse_lyso_withtree <- ordered(factor_exp5_lyso_withtree, levels = c("HO_LYS","HO_VEH","NO_LYS", "NO_VEH"))

###Plot###
dpi=600
tiff("Figure S7B UniFrac.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.5, 0.8), c(-0.3, 0.4), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(ord_lyso$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[factor_exp5_lys_withtree], lwd = 1)
points(ord_lyso$vectors[,1:2], pch = 21, cex = 1.3, col = pca_colors[factor_exp5_veh_withtree], lwd = 1)
ordiellipse(ord_lyso$vectors[,1:2], factor_exp5, kind = "se", conf = 0.95,col = pca_colors)
ordispider(ord$vectors[,1:2], factor_exp5, label = TRUE)
dev.off()


########################## Figure S7C ##########################

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



########################## Figure S7D-S7G ##########################

diff = phyloseq_to_deseq2(exp5_lyso_gen_prev, ~ SampleType)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")

###Pull out pairwise comparisons###
deseq2resultslyso2 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","HO_VEH","HO_LYS"))
deseq2resultslyso3 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","NO_LYS","HO_LYS"))
deseq2resultslyso4 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","NO_VEH","NO_LYS"))
deseq2resultslyso6 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","NO_VEH","HO_VEH"))

###Creates a data frame with the significantly differentially abundant genera###
diffabund_deseq2 <- function(deseq2results){
  deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
  deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(deseq2results), ], "matrix"))
  sigtab = deseq2results[(deseq2results$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(sigtab), ], "matrix"))
  sigtabgenus = subset(sigtab, !is.na(Genus))
  print(sigtabgenus)
}

deseq2lyso2 <- diffabund_deseq2(deseq2resultslyso2)
#Save .csv
deseq2lyso3 <- diffabund_deseq2(deseq2resultslyso3)
#Save .csv
deseq2lyso4 <- diffabund_deseq2(deseq2resultslyso4)
#Save .csv
deseq2lyso6 <- diffabund_deseq2(deseq2resultslyso6)
#save .csv

###Plots results###
deseq2plot <- function(sigtabgenus){
  diffabund <- ggplot(sigtabgenus,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
    geom_col() +
    coord_flip()
  diffabund2 <- diffabund + theme_bw()
  diffabund2
}

deseq2plot(deseq2lyso2)
#Save Figure S7F

deseq2plot(deseq2lyso3)
#Save Figure S7G

deseq2plot(deseq2lyso4)
#Save Figure S7D

deseq2plot(deseq2lyso6)
#Save Figure S7E


########################## Figure S7H ##########################

###Prepare data for barplot###
exp5_lyso_gen <- tax_glom(exp5_lyso, taxrank = 'Genus')
exp5_lyso_gen_prev <- filter_taxa(exp5_lyso_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)
exp5_lyso_gen_rel <- transform_sample_counts(exp5_lyso_gen_prev, function(x) x / sum(x) )
exp5_lyso_gen_df <- psmelt(exp5_lyso_gen_rel)

###Barplot with legend###
barplot_gen_lyso <- ggplot(exp5_lyso_gen_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~SampleType, scales = "free_x") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_lyso
#Save

###Barplot without legend###
barplot_gen_lyso_nolegend <- ggplot(exp5_lyso_gen_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~SampleType, scales = "free_x") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_lyso_nolegend
#Save
