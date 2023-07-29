library(tidyverse)
library(vegan)
library(phyloseq)
library(data.table)
library(mia)
library(scater)

########################## Figure S2A ##########################

#Load hyperoxia unfiltered OTU table with blanks
pankaj2_unfilt <- pankaj2 %>%
  tibble::column_to_rownames("Name")

#Load unfiltered taxa table with blanks
pankaj2.taxa.unfilt <- pankaj2.taxa.unfilt %>%
  tibble::column_to_rownames("Name")
pankaj2.taxa.unfilt <- as.matrix(pankaj2.taxa.unfilt)

#Load metadata with blanks

pankaj2.sort.unfilt <- pankaj2[,order(colnames(pankaj2))]
md.pankaj2.unfilt <- roughmetadata_pankaj2[which(unlist(roughmetadata_pankaj2$Name) %in% colnames(pankaj2.sort)),]

md.pankaj2.sort.unfilt <- md.pankaj2[order(unlist(md.pankaj2.unfilt$Name)),]
md.pankaj2.sort.unfilt <- md.pankaj2.sort.unfilt %>%
  tibble::column_to_rownames("Name")

summary(colSums(pankaj2.sort[,-1]))

OTU = otu_table(pankaj2.sort, taxa_are_rows = TRUE)
TAX = tax_table(pankaj2.taxa)
samples = sample_data(md.pankaj2.sort)

exp5_pankaj2 <- phyloseq(OTU, TAX, samples)

###Subset to just the blanks###
exp5_pankaj2_gen <- tax_glom(exp5_pankaj2, taxrank = 'Genus')
exp5_pankaj2_blank <- subset_samples(exp5_pankaj2_gen, Treatment == "Blank")

###Get top 50 OTUs from the blanks###
top20_pankaj2_blank_list <- names(sort(taxa_sums(exp5_pankaj2_blank), decreasing=TRUE)[1:50])
top20_pankaj2_blank_rel <- transform_sample_counts(exp5_pankaj2_blank, function(x) x / sum(x) )
top20_pankaj2_blank_df <- psmelt(top20_pankaj2_blank_rel)
top20_pankaj2_blank_df[!(top20_pankaj2_blank_df$OTU %in% top20_pankaj2_blank_list),]$Genus <- 'Other'

###Get top 50 OTUS from the filtered OTU table###
top20_pankaj2_gen_prev_list <- names(sort(taxa_sums(exp5_pankaj2_gen_prev), decreasing=TRUE)[1:50])
top20_pankaj2_gen_prev_rel <- transform_sample_counts(exp5_pankaj2_gen_prev, function(x) x / sum(x) )
top20_pankaj2_gen_prev_df <- psmelt(top20_pankaj2_gen_prev_rel)
top20_pankaj2_gen_prev_df[!(top20_pankaj2_gen_prev_df$OTU %in% top20_pankaj2_gen_prev_list),]$Genus <- 'Other'

#Merge!
blank_vs_samples <- rbind.fill(top20_pankaj2_blank_df, top20_pankaj2_gen_prev_df)

###Barplot with cleaned samples and blanks###

barplot_gen_pankaj2 <- ggplot(blank_vs_samples, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~Oxygen, scales = "free") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_pankaj2
#save


########################## Figure S2B ##########################

pankaj2.sort.flip.unfilt <- t(pankaj2.sort.unfilt)
simp_unfilt <- vegan::diversity(pankaj2.sort.flip.unfilt, "simpson")
shan_unfilt <- vegan::diversity(pankaj2.sort.flip.unfilt, "shannon")
richness_unfilt <- specnumber(pankaj2.sort.flip.unfilt)
pielou_unfilt <- shan/log(richness)
indices_unfilt <- as.data.frame(cbind(simp_unfilt, shan_unfilt,richness_unfilt,pielou_unfilt))
colnames(indices_unfilt) <- c("Simpson","Shannon","Species Richness","Pielou")

wilcox_alpha_pankaj2_unfilt <- t(sapply(indices_unfilt, function(x) unlist(kruskal.test(x~md.pankaj2.sort$Treatment)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_pankaj2_unfilt

richness_est_pankaj2_unfilt <- indices %>%
  mutate(
    Treatment = md.pankaj2.sort.unfilt$Treatment
  )

ggboxplot(richness_est_pankaj2_unfilt, x = "Treatment", y = "Species Richness",add = "jitter")
#Save


########################## Figure S2C ##########################

pankaj2.sort.flip <- t(pankaj2.sort)
simp <- vegan::diversity(pankaj2.sort.flip, "simpson")
shan <- vegan::diversity(pankaj2.sort.flip, "shannon")
richness <- specnumber(pankaj2.sort.flip)
pielou <- shan/log(richness)
indices <- as.data.frame(cbind(simp, shan,richness,pielou))
colnames(indices) <- c("Simpson","Shannon","Species Richness","Pielou")

wilcox_alpha_pankaj2 <- t(sapply(indices, function(x) unlist(kruskal.test(x~md.pankaj2.sort$Batch)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_pankaj2

richness_est_pankaj2 <- indices %>%
  mutate(
    Batch = md.pankaj2.sort$Batch
  )

#Compare Richness of Unfiltered vs Filtered samples
filt_vs_unfilt <- as.data.frame(cbind(richness_est_pankaj2[c(1:22),3],richness_est_pankaj2[,3]))
colnames(filt_vs_unfilt) <- c("Unfiltered","Filtered")
filt_unfilt_long <- pivot_longer(filt_vs_unfilt, cols = c("Unfiltered","Filtered"), names_to = "microDecon", values_to = "Richness")

ggboxplot(filt_unfilt_long, x = "microDecon", y = "Richness",add = "jitter")
#Save


########################## Figure S2D ##########################

#Load lysozyme unfiltered OTU table with blanks
lyso_unfilt <- lyso %>%
  tibble::column_to_rownames("Name")

#Load unfiltered taxa table with blanks
lyso.taxa.unfilt <- lyso.taxa.unfilt %>%
  tibble::column_to_rownames("Name")
lyso.taxa.unfilt <- as.matrix(lyso.taxa.unfilt)

#Load metadata with blanks

lyso.sort.unfilt <- lyso[,order(colnames(lyso))]
md.lyso.unfilt <- roughmetadata_lyso[which(unlist(roughmetadata_lyso$Name) %in% colnames(lyso.sort)),]

md.lyso.sort.unfilt <- md.lyso[order(unlist(md.lyso.unfilt$Name)),]
md.lyso.sort.unfilt <- md.lyso.sort.unfilt %>%
  tibble::column_to_rownames("Name")

summary(colSums(lyso.sort[,-1]))

OTU.lyso.unfilt = otu_table(lyso.sort.unfilt, taxa_are_rows = TRUE)
TAX.lyso.unfilt = tax_table(lyso.taxa.unfilt)
samples.lyso.unfilt = sample_data(md.lyso.sort.unfilt)

exp5_lyso_unfilt <- phyloseq(OTU, TAX, samples)


###Subset to just the blanks###
exp5_lyso_gen_unfilt <- tax_glom(exp5_lyso_unfilt, taxrank = 'Genus')
exp5_lyso_blank <- subset_samples(exp5_lyso_gen_unfilt, Treatment == "Blank")

###Get top 50 OTUs from the blanks###
top20_lyso_blank_list <- names(sort(taxa_sums(exp5_lyso_blank), decreasing=TRUE)[1:50])
top20_lyso_blank_rel <- transform_sample_counts(exp5_lyso_blank, function(x) x / sum(x) )
top20_lyso_blank_df <- psmelt(top20_lyso_blank_rel)
top20_lyso_blank_df[!(top20_lyso_blank_df$OTU %in% top20_lyso_blank_list),]$Genus <- 'Other'

###Get top 50 OTUS from the filtered OTU table###
top20_lyso_gen_prev_list <- names(sort(taxa_sums(exp5_lyso_gen_prev), decreasing=TRUE)[1:50])
top20_lyso_gen_prev_rel <- transform_sample_counts(exp5_lyso_gen_prev, function(x) x / sum(x) )
top20_lyso_gen_prev_df <- psmelt(top20_lyso_gen_prev_rel)
top20_lyso_gen_prev_df[!(top20_lyso_gen_prev_df$OTU %in% top20_lyso_gen_prev_list),]$Genus <- 'Other'

#Merge!
blank_vs_samples_lyso <- rbind.fill(top20_lyso_blank_df, top20_lyso_gen_prev_df)

###Barplot with cleaned samples and blanks###

barplot_gen_lyso_blank <- ggplot(blank_vs_samples_lyso, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~SampleType, scales = "free") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_lyso_blank
#save

########################## Figure S2E ##########################

lyso.sort.flip.unfilt <- t(lyso.sort.unfilt)
simp_lyso_unfilt <- vegan::diversity(lyso.sort.flip.unfilt, "simpson")
shan_lyso_unfilt <- vegan::diversity(lyso.sort.flip.unfilt, "shannon")
richness_lyso_unfilt <- specnumber(lyso.sort.flip.unfilt)
pielou_lyso_unfilt <- shan/log(richness)
indices_lyso_unfilt <- as.data.frame(cbind(simp_lyso_unfilt, shan_lyso_unfilt,richness_lyso_unfilt,pielou_lyso_unfilt))
colnames(indices_lyso_unfilt) <- c("Simpson","Shannon","Species Richness","Pielou")

wilcox_alpha_lyso_unfilt <- t(sapply(indices_unfilt, function(x) unlist(kruskal.test(x~md.lyso.sort$Treatment)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_lyso_unfilt

richness_est_lyso_unfilt <- indices_lyso %>%
  mutate(
    Treatment = md.lyso.sort.unfilt$Treatment
  )

ggboxplot(richness_est_lyso_unfilt, x = "Treatment", y = "Species Richness",add = "jitter")
#Save


########################## Figure S2F ##########################

###This part is the filtered data###
lyso.sort.flip <- t(lyso.sort)
simp_lyso <- vegan::diversity(lyso.sort.flip, "simpson")
shan_lyso <- vegan::diversity(lyso.sort.flip, "shannon")
richness_lyso <- specnumber(lyso.sort.flip)
pielou_lyso <- shan/log(richness)
indices_lyso <- as.data.frame(cbind(simp_lyso, shan_lyso,richness_lyso,pielou_lyso))
colnames(indices_lyso) <- c("Simpson","Shannon","Species Richness","Pielou")

wilcox_alpha_lyso <- t(sapply(indices_lyso, function(x) unlist(kruskal.test(x~md.lyso.sort$Batch)[c("estimate","p.value","statistic","conf.int")])))
wilcox_alpha_lyso

richness_est_lyso <- indices %>%
  mutate(
    Batch = md.lyso.sort$Batch
  )

#Compare Richness of Unfiltered vs Filtered samples
filt_vs_unfilt_lyso <- as.data.frame(cbind(richness_est_lyso_unfilt[c(1:22),3],richness_est_lyso[,3]))
colnames(filt_vs_unfilt_lyso) <- c("Unfiltered","Filtered")
filt_unfilt_long_lyso <- pivot_longer(filt_vs_unfilt_lyso, cols = c("Unfiltered","Filtered"), names_to = "microDecon", values_to = "Richness")

ggboxplot(filt_unfilt_long, x = "microDecon", y = "Richness",add = "jitter")
#Save

########################## Figure S2G ##########################

exp5_pankaj2_unfilt_noblank <- subset_samples(exp5_pankaj2_gen, Batch != "Blank")
pankaj2_unfilt <- makeTreeSummarizedExperimentFromPhyloseq(exp5_pankaj2_unfilt_noblank)

prevalence.frequency <- getPrevalence(pankaj2_unfilt,
                                      rank = "Genus",
                                      detection = 0,
                                      sort = TRUE,
                                      as_relative = TRUE)

rowData(pankaj2_unfilt)$prevalence <- prevalence.frequency
#Prevalence plot without legend
plotRowData(pankaj2_unfilt, "prevalence", colour_by = "Genus", add_legend = FALSE)
#Save

#Prevalence plot with legend
plotRowData(pankaj2_unfilt, "prevalence", colour_by = "Genus", add_legend = TRUE)
#Save
