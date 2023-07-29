library(tidyverse)
library(readxl)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(factoextra)
library(microDecon)

###Load Unfiltered Data with Blanks###
#Load OTUs
#Load Taxa
#Load Metadata

#Note: OTUs ordered according to experiment batch, blanks are the first two samples listed. Taxa laid out as a single column with semicolons between taxonomic levels

#Set everything up
bact <- bact %>%
  tibble::column_to_rownames("Name")
bact[is.na(bact)] <- 0

bact.taxa <- bact.taxa %>%
  tibble::column_to_rownames("Name")
bact.taxa <- as.matrix(bact.taxa)

md.bact <- roughmetadata_bact[which(unlist(roughmetadata_bact$Name) %in% colnames(bact)),]
summary(colSums(bact[,-1]))

md.bact <- md.bact %>%
  tibble::column_to_rownames("Name")

OTU_rough = otu_table(bact, taxa_are_rows = TRUE)
TAX_rough = tax_table(bact.taxa)
samples_rough = sample_data(md.bact)

#Convert to phyloseq object
exp1_bact_rough <- phyloseq(OTU_rough, TAX_rough, samples_rough)

#Extract components from phyloseq object and set everything up in the proper format for microDecon
exp1_bact_otu <- as.data.frame(exp1_bact_rough@otu_table)
exp1_bact_tax <- as.data.frame(exp1_bact_rough@tax_table)
exp1_bact_meta <- as.data.frame(exp1_bact_rough@sam_data)

exp1_bact_otu <- tibble::rownames_to_column(exp1_bact_otu, "Name")
exp1_bact_microdecon_format <- data.frame(c(exp1_bact_otu,exp1_bact_tax))
#Save .csv, this OTU table can be used later to compare with the decontaminated table

#Run microDecon and save output
set.seed(1312)
decon_try1 <- decon(data = exp1_bact_microdecon_format,numb.blanks=2,numb.ind=c(14,8),taxa=T, runs=2)

results_try1 <- decon_try1$decon.table
#Save as .csv

removed_try1 <- decon_try1$reads.removed
#Save as .csv

mean_try1 <- decon_try1$mean.per.group
#Save as .csv

sum_try1 <- decon_try1$sum.per.group
#Save as .csv

asvs_try1 <- decon_try1$OTUs.removed
#Save as .csv
