library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(data.table)
library(RColorBrewer)
library(Maaslin2)
library(pulsar)
library(huge)
library(MASS)
library(VGAM)
library(Matrix)
library(testthat)
library(parallel)
library(boot)
library(SpiecEasi)
library(igraph)
library(RColorBrewer)
library(randomForest)
library(xgboost)
library(e1071)
library(caret)
library(glmnet)
library(MASS)
library(pROC)
library(stringr)
library(MLeval)

########Barplot########

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
ggsave("2023_01_19_exp5_1_barplot_legend_2.png", width = 15, height = 10, device='png', dpi=600)

###Barplot without legend###

barplot_gen_lyso_nolegend <- ggplot(exp5_lyso_gen_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~SampleType, scales = "free_x") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_lyso_nolegend
ggsave("2023_01_19_exp5_1_barplot_nolegend_2.png", width = 15, height = 10, device='png', dpi=600)


########Genus Level DESeq2 Analysis########

exp5_lyso_gen <- tax_glom(exp5_lyso, taxrank = 'Genus')
exp5_lyso_gen_prev <- filter_taxa(exp5_lyso_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)

diff = phyloseq_to_deseq2(exp5_lyso_gen_prev, ~ SampleType)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diff), 1, gm_mean)
diff = estimateSizeFactors(diff, geoMeans = geoMeans)
diff = DESeq(diff, fitType="local")

###Separate out pairwise comparisons###
deseq2resultslyso1 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_VEH","HO_LYS"))
deseq2resultslyso2 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","HO_VEH","HO_LYS"))
deseq2resultslyso3 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_LYS","HO_LYS"))
deseq2resultslyso4 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_VEH","RA_LYS"))
deseq2resultslyso5 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","HO_VEH","RA_LYS"))
deseq2resultslyso6 = results(diff, pAdjustMethod = "fdr",contrast = c("SampleType","RA_VEH","HO_VEH"))

###Create dataframe of significant differences###
diffabund_deseq2_gen <- function(deseq2results){
  deseq2results = deseq2results[order(deseq2results$padj, na.last=NA), ]
  deseq2results_df = cbind(as(deseq2results, "data.frame"), as(tax_table(exp5_lyso_gen_prev)[rownames(deseq2results), ], "matrix"))
  sigtab = deseq2results[(deseq2results$padj < 0.05), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(exp5_lyso_gen)[rownames(sigtab), ], "matrix"))
  sigtabgenus = subset(sigtab, !is.na(Genus))
  return(sigtabgenus)
}

deseq2lyso1 <- diffabund_deseq2_gen(deseq2resultslyso1)
write.csv(deseq2lyso1, "2022_01_19_exp5_lyso_diffabund_deseq_raveh_holys_stats.csv")
deseq2lyso2 <- diffabund_deseq2_gen(deseq2resultslyso2)
write.csv(deseq2lyso2, "2022_01_19_exp5_lyso_diffabund_deseq_hoveh_holys_stats.csv")
deseq2lyso3 <- diffabund_deseq2_gen(deseq2resultslyso3)
write.csv(deseq2lyso3, "2022_01_19_exp5_lyso_diffabund_deseq_ralys_holys_stats.csv")
deseq2lyso4 <- diffabund_deseq2_gen(deseq2resultslyso4)
write.csv(deseq2lyso4, "2022_01_19_exp5_lyso_diffabund_deseq_raveh_ralys_stats.csv")
deseq2lyso5 <- diffabund_deseq2_gen(deseq2resultslyso5)
write.csv(deseq2lyso5, "2022_01_19_exp5_lyso_diffabund_deseq_hoveh_ralys_stats.csv")
deseq2lyso6 <- diffabund_deseq2_gen(deseq2resultslyso6)
write.csv(deseq2lyso6, "2022_01_19_exp5_lyso_diffabund_deseq_raveh_hoveh_stats.csv")

###Plot results!###
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

########SPIEC-EASI########

#Filter to only OTUs found in 50% of samples and not listed as NA at the Genus level
exp5_lyso_prev <- filter_taxa(exp5_lyso , function(x) sum(x >= 1) > (0.5*length(x)), TRUE)
exp5_lyso_prev_narm = subset_taxa(exp5_lyso_prev, Genus!="NA")

#Pull out the OTU, taxonomy, and metadata tables as data frames
exp5_lyso_narm_otu <- as.data.frame(exp5_lyso_prev_narm@otu_table)
exp5_lyso_narm_tax <- as.data.frame(exp5_lyso_prev_narm@tax_table)
exp5_lyso_narm_meta <- as.data.frame(exp5_lyso_prev_narm@sam_data)

exp5_lyso_narm_otu_ho_lys <- exp5_lyso_narm_otu[,which(exp5_lyso_narm_meta$SampleType == "HO_LYS")]
exp5_lyso_narm_otu_ho_veh <- exp5_lyso_narm_otu[,which(exp5_lyso_narm_meta$SampleType == "HO_VEH")]
exp5_lyso_narm_otu_ra_lys <- exp5_lyso_narm_otu[,which(exp5_lyso_narm_meta$SampleType == "RA_LYS")]
exp5_lyso_narm_otu_ra_veh <- exp5_lyso_narm_otu[,which(exp5_lyso_narm_meta$SampleType == "RA_VEH")]

#ho_lys
exp5_lyso_narm_tax<-as.matrix(exp5_lyso_narm_tax)
OTU_lyso_narm_ho_lys = otu_table(exp5_lyso_narm_otu_ho_lys, taxa_are_rows = TRUE)
TAX_lyso_narm_ho_lys= tax_table(exp5_lyso_narm_tax)
samples_lyso_narm_ho_lys= sample_data(exp5_lyso_narm_meta)

exp5_lyso_narm_ho_lys <- phyloseq(OTU_lyso_narm_ho_lys, TAX_lyso_narm_ho_lys, samples_lyso_narm_ho_lys)

#ho_veh
exp5_lyso_narm_tax<-as.matrix(exp5_lyso_narm_tax)
OTU_lyso_narm_ho_veh = otu_table(exp5_lyso_narm_otu_ho_veh, taxa_are_rows = TRUE)
TAX_lyso_narm_ho_veh= tax_table(exp5_lyso_narm_tax)
samples_lyso_narm_ho_veh= sample_data(exp5_lyso_narm_meta)

exp5_lyso_narm_ho_veh <- phyloseq(OTU_lyso_narm_ho_veh, TAX_lyso_narm_ho_veh, samples_lyso_narm_ho_veh)

#RA_LYS
exp5_lyso_narm_tax<- as.matrix(exp5_lyso_narm_tax)
OTU_lyso_narm_ra_lys = otu_table(exp5_lyso_narm_otu_ra_lys, taxa_are_rows = TRUE)
TAX_lyso_narm_ra_lys= tax_table(exp5_lyso_narm_tax)
samples_lyso_narm_ra_lys= sample_data(exp5_lyso_narm_meta)

exp5_lyso_narm_ra_lys <- phyloseq(OTU_lyso_narm_ra_lys, TAX_lyso_narm_ra_lys, samples_lyso_narm_ra_lys)

#RA_VEH
exp5_lyso_narm_tax<- as.matrix(exp5_lyso_narm_tax)
OTU_lyso_narm_ra_veh = otu_table(exp5_lyso_narm_otu_ra_veh, taxa_are_rows = TRUE)
TAX_lyso_narm_ra_veh= tax_table(exp5_lyso_narm_tax)
samples_lyso_narm_ra_veh= sample_data(exp5_lyso_narm_meta)

exp5_lyso_narm_ra_veh <- phyloseq(OTU_lyso_narm_ra_veh, TAX_lyso_narm_ra_veh, samples_lyso_narm_ra_veh)

#Run SPIEC-EASI, set threshold at 0.1 as recommended                          
se.exp5.lyso.ho.lys <- spiec.easi(exp5_lyso_narm_ho_lys, method='mb', nlambda=99,
                                  lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp5.lyso.ho.veh <- spiec.easi(exp5_lyso_narm_ho_veh, method='mb', nlambda=99,
                                  lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp5.lyso.ra.lys <- spiec.easi(exp5_lyso_narm_ra_lys, method='mb', nlambda=99,
                                  lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp5.lyso.ra.veh <- spiec.easi(exp5_lyso_narm_ra_veh, method='mb', nlambda=99,
                                  lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

#HO_LYS Plot
spiec.graph.lyso.ho.lys=adj2igraph(getRefit(se.exp5.lyso.ho.lys), vertex.attr=list(name=taxa_names(exp5_lyso_narm_ho_lys)))

narm_spiec_plot_lyso_ho_lys <- plot_network(spiec.graph.lyso.ho.lys, exp5_lyso_narm_ho_lys, type='taxa', color="Genus")
narm_spiec_plot_lyso_ho_lys + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_lyso_spieceasi_ho_lys.png", width = 10, height = 5, device='png', dpi=600)

#How many OTUs are in the network
nodes_ho_lys <- gorder(spiec.graph.lyso.ho.lys)
#How many connections between OTUs are in the network
edges_ho_lys <- gsize(spiec.graph.lyso.ho.lys)
#Number of shortest paths through node, measure of centrality
betweenness_ho_lys <- betweenness(spiec.graph.lyso.ho.lys, normalized = TRUE)
#How many edges per node
degree_ho_lys <- as.data.frame(degree(spiec.graph.lyso.ho.lys))
degree_distribution_ho_lys <- degree_distribution(spiec.graph.lyso.ho.lys)

#HO_VEH Plot
spiec.graph.lyso.ho.veh=adj2igraph(getRefit(se.exp5.lyso.ho.veh), vertex.attr=list(name=taxa_names(exp5_lyso_narm_ho_veh)))

narm_spiec_plot_lyso_ho_veh <- plot_network(spiec.graph.lyso.ho.veh, exp5_lyso_narm_ho_veh, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_lyso_ho_veh + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_lyso_spieceasi_ho_veh.png", width = 10, height = 5, device='png', dpi=600)

nodes_ho_veh <- gorder(spiec.graph.lyso.ho.veh)
edges_ho_veh <- gsize(spiec.graph.lyso.ho.veh)
betweenness_ho_veh <- as.list(betweenness(spiec.graph.lyso.ho.veh, normalized = TRUE))
degree_ho_veh <- as.data.frame(degree(spiec.graph.lyso.ho.veh))


#RA_LYS Plot
spiec.graph.lyso.ra.lys=adj2igraph(getRefit(se.exp5.lyso.ra.lys), vertex.attr=list(name=taxa_names(exp5_lyso_narm_ra_lys)))

narm_spiec_plot_lyso_ra_lys <- plot_network(spiec.graph.lyso.ra.lys, exp5_lyso_narm_ra_lys, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_lyso_ra_lys + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_lyso_spieceasi_ra_lys.png", width = 10, height = 5, device='png', dpi=600)

nodes_ra_lys <- gorder(spiec.graph.lyso.ra.lys)
edges_ra_lys <- gsize(spiec.graph.lyso.ra.lys)
betweenness_ra_lys <- as.list(betweenness(spiec.graph.lyso.ra.lys, normalized = TRUE))
degree_ra_lys <- as.data.frame(degree(spiec.graph.lyso.ra.lys))

#RA_VEH Plot
spiec.graph.lyso.ra.veh=adj2igraph(getRefit(se.exp5.lyso.ra.veh), vertex.attr=list(name=taxa_names(exp5_lyso_narm_ra_veh)))

narm_spiec_plot_lyso_ra_veh <- plot_network(spiec.graph.lyso.ra.veh, exp5_lyso_narm_ra_veh, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_lyso_ra_veh + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_lyso_spieceasi_ra_veh.png", width = 10, height = 5, device='png', dpi=600)


nodes_ra_veh <- gorder(spiec.graph.lyso.ra.veh)
edges_ra_veh <- gsize(spiec.graph.lyso.ra.veh)
betweenness_ra_veh <- as.list(betweenness(spiec.graph.lyso.ra.veh, normalized = TRUE))
degree_ra_veh <- as.data.frame(degree(spiec.graph.lyso.ra.veh))

###Degree plot###
degree_df <- data.frame(degree_ho_lys,degree_ho_veh,degree_ra_lys,degree_ra_veh)
colnames(degree_df) <- c("HO_LYS","HO_VEH","RA_LYS","RA_VEH")

count_ho_lys <- as.data.frame(table(degree_df$HO_LYS))
colnames(count_ho_lys) <- c("Degree","HO_LYS")

count_ho_veh <- as.data.frame(table(degree_df$HO_VEH))
colnames(count_ho_veh) <- c("Degree","HO_VEH")
#Adds an empty row so column lengths match
count_ho_veh[nrow(count_ho_veh) + 1,] <- c("8",0)

count_ra_lys <- as.data.frame(table(degree_df$RA_LYS))
colnames(count_ra_lys) <- c("Degree","RA_LYS")
count_ra_lys[nrow(count_ra_lys) + 1,] <- c("8",0)

count_ra_veh <- as.data.frame(table(degree_df$RA_VEH))
colnames(count_ra_veh) <- c("Degree","RA_VEH")

#Plots how many nodes have a given degree for each group
degree_freq <- data.frame(count_ho_lys, count_ho_veh$HO_VEH, count_ra_lys$RA_LYS, count_ra_veh$RA_VEH)
colnames(degree_freq) <- c("Degree","HO_LYS","HO_VEH","RA_LYS","RA_VEH")
degree_freq$HO_VEH <- as.integer(degree_freq$HO_VEH)
degree_freq$RA_LYS <- as.integer(degree_freq$RA_LYS)

degree_long <- degree_freq %>%
  pivot_longer(c("HO_LYS","HO_VEH","RA_LYS","RA_VEH"), names_to = "SampleType", values_to = "Frequency")


ggplot(degree_long,aes(x = Degree,y = Frequency, group = SampleType))+
  geom_line()+
  theme(color = SampleType)

ggplot(degree_long,aes(x=Degree,y=Frequency,group=SampleType))+
  geom_line(aes(color=SampleType)) +
  geom_point(aes(color=SampleType)) +
  theme_bw()
ggsave("2023_01_23_exp5_lyso_spieceasi_degree.png", width = 5, height = 5, device='png', dpi=600)

###Edges Plot###
edges_df <- data.frame(edges_ho_lys,edges_ho_veh,edges_ra_lys,edges_ra_veh)
colnames(edges_df) <- c("HO_LYS","HO_VEH","RA_LYS","RA_VEH")
edges_long <- edges_df %>%
  pivot_longer(c("HO_LYS","HO_VEH","RA_LYS","RA_VEH"), names_to = "SampleType", values_to = "Edges")

#Plots the number of edges in each group
ggplot(edges_long,aes(x= SampleType, y = Edges))+
  geom_bar(stat="identity", aes(fill=SampleType))+
  theme_bw()
ggsave("2023_01_23_exp5_lyso_spieceasi_edges.png", width = 5, height = 5, device='png', dpi=600)  


########Feature Selection (Random Forest)########

exp5_lyso_prev_noNA = subset_taxa(exp5_lyso_gen_prev, Genus!="NA")
exp5_lyso_prev_rel <- transform_sample_counts(exp5_lyso_prev_noNA, function(x) x / sum(x) )

#Prep OTU table
predictors <- as.data.frame(exp5_lyso_prev_rel@otu_table)
taxlabels <- as.data.frame(exp5_lyso_prev_rel@tax_table)
predictors <- t(predictors)
colnames(predictors) <- taxlabels$Genus
#Get relevant metadata
response <- as.factor(sample_data(exp5_lyso_prev_rel)$SampleType)

#Make dataframe combining OTUs and Metadata
rf.data <- data.frame(response, predictors)

set.seed(1312)
lyso_idx = createDataPartition(rf.data$response, p = 0.75, list = FALSE)
lyso_train = rf.data[lyso_idx, ]
lyso_test = rf.data[-lyso_idx, ]

#Leave One Out Cross-Validation
lyso_cv <- trainControl(method='LOOCV', savePredictions = TRUE)

#Tell model to weight all groups evenly
model_weights <- ifelse(lyso_train$response == "HO_LYS",  (1/table(lyso_train$response)[1]) * 0.25,
                        ifelse(lyso_train$response == "HO_VEH",  (1/table(lyso_train$response)[2]) * 0.25,
                               ifelse(lyso_train$response == "RA_LYS",  (1/table(lyso_train$response)[3]) * 0.25,
                                      (1/table(lyso_train$response)[4]) * 0.25))) 

###Set parameter tuning grid, testing a few parameter combinations###
tgrid <- expand.grid(
  .mtry = c(2,16,32,63),
  .splitrule = c("gini","extratrees"),
  .min.node.size = c(1,3,5,10,18)
)

###Run random forest###
set.seed(1312)
rf_fit <- train(as.factor(response) ~ ., 
                data = lyso_train, 
                method = "ranger",
                tuneGrid = tgrid,
                num.trees = 1000,
                trControl = lyso_cv,
                weights = model_weights,
                verbose = FALSE,
                verbosity = 0,
                importance = "permutation"
)

rf_fit

###New grid, using the best hyperparameter combination###
best_grid <- expand.grid(
  .mtry = 63,
  .splitrule = "extratrees",
  .min.node.size = 5
)

###Final random forest model###
set.seed(1312)
rf_fit2 <- train(as.factor(response) ~ ., 
                 data = lyso_train, 
                 method = "ranger",
                 tuneGrid = best_grid,
                 num.trees = 1000,
                 trControl = lyso_cv,
                 weights = model_weights,
                 metric = "Sens",
                 verbose = FALSE,
                 verbosity = 0,
                 importance = "permutation"
)

rf_fit2

###Plot ROC###
res_rf <- evalm(rf_fit2)
ggsave("lyso_randomforest_roc.png", width = 5, height = 5, device='png', dpi=600)

###Check against test data###
testclass <- predict(rf_fit2, newdata = lyso_test)
cfMatrix <- confusionMatrix(data = testclass, lyso_test$response)

###Plot feature importance###
features <- varImp(rf_fit2)
features <- features$importance
features <- rownames_to_column(features)
features <- features[order(-features$Overall),]
features_top50 <- as.data.frame(features[1:50,])

ggplot2::ggplot(features_top50, aes(x=rowname, y=Overall)) +
  geom_point( color="#FB0207", size=4, alpha=1.0)+
  geom_segment( aes(x=rowname, xend=rowname, y=0, yend=Overall), 
                color='#EF94A2') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light() +
  coord_flip() 

ggsave("lyso_randomforest.png", width = 5, height = 10, device='png', dpi=600)
