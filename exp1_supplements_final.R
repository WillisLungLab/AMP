###Supplementary figures###

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
exp5_pankaj_gen <- tax_glom(exp5_pankaj, taxrank = 'Genus')
exp5_pankaj_gen_prev <- filter_taxa(exp5_pankaj_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)
exp5_pankaj_gen_rel <- transform_sample_counts(exp5_pankaj_gen_prev, function(x) x / sum(x) )
exp5_pankaj_gen_df <- psmelt(exp5_pankaj_gen_rel)

###Barplot of top 20 combined genera, and an entry for all others###

barplot_gen_pankaj <- ggplot(exp5_pankaj_gen_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~SampleType, scales = "free_x") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_pankaj
ggsave("2023_01_19_exp5_1_barplot_legend_2.png", width = 15, height = 10, device='png', dpi=600)


########Differential Abundance Barplots########

###DESeq2 Supplements###
diffabund_deseq2 <- ggplot(sigtab_otu,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
  geom_col() +
  coord_flip()+
  scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
diffabund_deseq2 + theme_bw()

###MaAsLin2 Supplements###
sigtabmaaslin2 <- read.delim2("exp5_pankaj_output/significant_results.tsv")
diffabund_maaslin2 <- ggplot(sigtabmaaslin2,aes(x=feature,y=coef,fill=coef>0))+
  geom_col() + coord_flip()+
  scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
diffabund_maaslin2 + theme_bw()
ggsave("2023_01_19_exp5_1_diffabund_maaslin2.tiff", width = 5, height = 5, device='tiff', dpi=600)


########SpiecEasi########

#Filter to only OTUs found in 50% of samples and not listed as NA at the Genus level
exp5_pankaj_prev <- filter_taxa(exp5_pankaj , function(x) sum(x >= 1) > (0.50*length(x)), TRUE)
exp5_pankaj_prev_narm = subset_taxa(exp5_pankaj_prev, Genus!="NA")

###Pull out the OTU, taxonomy, and metadata tables as data frames###
exp5_pankaj_narm_otu <- as.data.frame(exp5_pankaj_prev_narm@otu_table)
exp5_pankaj_narm_tax <- as.data.frame(exp5_pankaj_prev_narm@tax_table)
exp5_pankaj_narm_meta <- as.data.frame(exp5_pankaj_prev_narm@sam_data)

###Separate out HO and RA###
exp5_pankaj_narm_otu_ho <- exp5_pankaj_narm_otu[,which(exp5_pankaj_narm_meta$SampleType == "HO")]
exp5_pankaj_narm_otu_ra <- exp5_pankaj_narm_otu[,which(exp5_pankaj_narm_meta$SampleType == "RA")]

###Prepare HO samples###
exp5_pankaj_narm_tax<-as.matrix(exp5_pankaj_narm_tax)
OTU_pankaj_narm_ho = otu_table(exp5_pankaj_narm_otu_ho, taxa_are_rows = TRUE)
TAX_pankaj_narm_ho= tax_table(exp5_pankaj_narm_tax)
samples_pankaj_narm_ho= sample_data(exp5_pankaj_narm_meta)

exp5_pankaj_narm_ho <- phyloseq(OTU_pankaj_narm_ho, TAX_pankaj_narm_ho, samples_pankaj_narm_ho)

###Prepare RA Samples###
exp5_pankaj_narm_tax<- as.matrix(exp5_pankaj_narm_tax)
OTU_pankaj_narm_ra = otu_table(exp5_pankaj_narm_otu_ra, taxa_are_rows = TRUE)
TAX_pankaj_narm_ra= tax_table(exp5_pankaj_narm_tax)
samples_pankaj_narm_ra= sample_data(exp5_pankaj_narm_meta)

exp5_pankaj_narm_ra <- phyloseq(OTU_pankaj_narm_ra, TAX_pankaj_narm_ra, samples_pankaj_narm_ra)

###Run SPIEC-EASI, set threshold at 0.1 as recommended###                        
se.exp5.pankaj.ho <- spiec.easi(exp5_pankaj_narm_ho, method='mb', nlambda=99,
                                lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

se.exp5.pankaj.ra <- spiec.easi(exp5_pankaj_narm_ra, method='mb', nlambda=99,
                                lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.1))

###HO Plot###
spiec.graph.pankaj.ho=adj2igraph(getRefit(se.exp5.pankaj.ho), vertex.attr=list(name=taxa_names(exp5_pankaj_narm_ho)))

narm_spiec_plot_pankaj_ho <- plot_network(spiec.graph.pankaj.ho, exp5_pankaj_narm_ho, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_pankaj_ho + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_1_spieceasi_ho.png", width = 10, height = 5, device='png', dpi=600)

narm_spiec_plot_pankaj_ho_nolegend <- plot_network(spiec.graph.pankaj.ho, exp5_pankaj_narm_ho, type='taxa', color="Genus", label = NULL)+
  theme(legend.position = "none")
narm_spiec_plot_pankaj_ho_nolegend + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_1_spieceasi_ho_nolegend.png", width = 10, height = 5, device='png', dpi=600)

#How many OTUs are present
nodes_ho <- gorder(spiec.graph.pankaj.ho)
#Number of connection between OTUs
edges_ho <- gsize(spiec.graph.pankaj.ho)
#Number of shortest paths through node, measure of centrality
betweenness_ho <- as.list(betweenness(spiec.graph.pankaj.ho, normalized = TRUE))
#Number of edges per node
degree_ho <- as.data.frame(degree(spiec.graph.pankaj.ho))

###RA Plot###
spiec.graph.pankaj.ra=adj2igraph(getRefit(se.exp5.pankaj.ra), vertex.attr=list(name=taxa_names(exp5_pankaj_narm_ra)))

narm_spiec_plot_pankaj_ra <- plot_network(spiec.graph.pankaj.ra, exp5_pankaj_narm_ra, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_pankaj_ra + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_1_spieceasi_ra.png", width = 10, height = 5, device='png', dpi=600)

narm_spiec_plot_pankaj_ra_nolegend <- plot_network(spiec.graph.pankaj.ra, exp5_pankaj_narm_ra, type='taxa', color="Genus", label = NULL)+
  theme(legend.position = "none")
narm_spiec_plot_pankaj_ra_nolegend + geom_point(shape = 1,size = 4,colour = "black")
ggsave("2023_01_21_exp5_1_spieceasi_ra_nolegend.png", width = 10, height = 5, device='png', dpi=600)

nodes_ra <- gorder(spiec.graph.pankaj.ra)
edges_ra <- gsize(spiec.graph.pankaj.ra)
betweenness_ra <- as.list(betweenness(spiec.graph.pankaj.ra, normalized = TRUE))
degree_ra <- as.data.frame(degree(spiec.graph.pankaj.ra))

###Degree plot###
degree_df <- data.frame(degree_ho,degree_ra)
colnames(degree_df) <- c("HO","RA")

count_ho <- as.data.frame(table(degree_df$HO))
colnames(count_ho) <- c("Degree","HO")

count_ra <- as.data.frame(table(degree_df$RA))
colnames(count_ra) <- c("Degree","RA")

degree_freq <- rbind.fill(count_ho[c("Degree", "HO")], count_ra[c("Degree", "RA")])

degree_freq <- degree_freq %>%
  group_by(Degree) %>%
  dplyr::summarise(across(c(HO,RA), ~sum(., na.rm = TRUE)))

degree_long <- degree_freq %>%
  pivot_longer(c("HO","RA"), names_to = "SampleType", values_to = "Frequency")

#Plot how many nodes have a certain degree for both groups
ggplot(degree_long,aes(x=Degree,y=Frequency,group=SampleType))+
  geom_line(aes(color=SampleType)) +
  geom_point(aes(color=SampleType)) +
  theme_bw()
ggsave("exp5_1_spieceasi_degreeplot_1_23_2023.png", width = 5, height = 5, device='png', dpi=600)

###Edges Plot###
edges_df <- data.frame(edges_ho,edges_ra)
colnames(edges_df) <- c("HO","RA")
edges_long <- edges_df %>%
  pivot_longer(c("HO","RA"), names_to = "SampleType", values_to = "Edges")

#Plot total edges in each group
ggplot(edges_long,aes(x= SampleType, y = Edges))+
  geom_bar(stat="identity", aes(fill=SampleType))+
  theme_bw()
ggsave("exp5_1_spieceasi_edgeplot_1_23_2023.png", width = 5, height = 5, device='png', dpi=600) 


########Random Forest########

exp5_combined_gen_prev_noNA = subset_taxa(exp5_pankaj_gen_prev, Genus!="NA")
exp5_combined_gen_prev_rel <- transform_sample_counts(exp5_combined_gen_prev_noNA, function(x) x / sum(x) )

###Prep OTU table###
predictors <- as.data.frame(exp5_combined_gen_prev_rel@otu_table)
taxlabels <- as.data.frame(exp5_combined_gen_prev_rel@tax_table)
predictors <- t(predictors)
colnames(predictors) <- taxlabels$Genus
#Get relevant metadata
response <- as.factor(sample_data(exp5_combined_gen_prev_rel)$SampleType)
#Surprise tool that will help us later
tax.labels <- as.data.frame(exp5_combined_gen_prev_rel@tax_table)

###Make dataframe combining OTUs and Metadata###
rf.data <- data.frame(response, predictors)

###Set random seed and partition data into training and testing###
set.seed(1312)
pankaj_idx = createDataPartition(rf.data$response, p = 0.75, list = FALSE)
pankaj_train = rf.data[pankaj_idx, ]
pankaj_test = rf.data[-pankaj_idx, ]

###Set cross-validation, in this case Leave One Out Cross-Validation###
pankaj_cv <- trainControl(method='LOOCV', classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

#Weight model so both groups are weighted equally
model_weights <- ifelse(pankaj_train$response == "HO",
                        (1/table(pankaj_train$response)[1]) * 0.5,
                        (1/table(pankaj_train$response)[2]) * 0.5
)

###Set parameter tuning grid, testing a few parameter combinations###
tgrid <- expand.grid(
  .mtry = c(2,16,32,63),
  .splitrule = c("gini","extratrees"),
  .min.node.size = c(1,3,5,10,18)
)

###Run random forest###
set.seed(1312)
rf_fit <- train(as.factor(response) ~ ., 
                data = pankaj_train, 
                method = "ranger",
                tuneGrid = tgrid,
                num.trees = 1000,
                trControl = pankaj_cv,
                weights = model_weights,
                metric = "ROC",
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
                data = pankaj_train, 
                method = "ranger",
                tuneGrid = best_grid,
                num.trees = 1000,
                trControl = pankaj_cv,
                weights = model_weights,
                metric = "Sens",
                verbose = FALSE,
                verbosity = 0,
                importance = "permutation"
)

rf_fit2

###Plot ROC###
res_rf <- evalm(rf_fit2)
ggsave("pankaj_randomforest_roc.png", width = 5, height = 5, device='png', dpi=600)

###Check against test data###
testclass <- predict(rf_fit2, newdata = pankaj_test)
cfMatrix <- confusionMatrix(data = testclass, pankaj_test$response)

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

ggsave("pankaj_randomforest.png", width = 5, height = 10, device='png', dpi=600)

