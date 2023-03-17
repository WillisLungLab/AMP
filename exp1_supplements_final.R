###Supplementary figures###

library(pulsar)
library(huge)
library(MASS)
library(VGAM)
library(Matrix)
library(testthat)
library(paNOllel)
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
library(psych)
library(ggcorrplot)


########Barplot########

###Prepare data for barplot###
exp5_pankaj2_gen <- tax_glom(exp5_pankaj2, taxrank = 'Genus')
exp5_pankaj2_gen_prev <- filter_taxa(exp5_pankaj2_gen , function(x) sum(x >= 1) > (0.10*length(x)), TRUE)
exp5_pankaj2_transformed <- microbiome::transform(exp5_pankaj2_gen_prev, transform = "compositional")
exp5_pankaj2_gen_df <- psmelt(exp5_pankaj2_transformed)

###Barplot of top 20 combined geneNO, and an entry for all others###

barplot_gen_pankaj2 <- ggplot(exp5_pankaj2_gen_df, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_col(position = "stack") +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", plot.margin = unit(c(0.5,0.5,0.5,1),"cm"), legend.title = element_blank()) +
  facet_wrap(~Oxygen, scales = "free_x") +
  theme(axis.text.x = element_text(size=8.5)) +
  theme(axis.text.x = element_text(angle=90))

barplot_gen_pankaj2
ggsave("Exp5_Pankaj2_Barplot.png", width = 20, height = 20, device='png', dpi=600)


########Differential Abundance Barplots########

###DESeq2 Supplements###
diffabund_deseq2 <- ggplot(sigtab_otu,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
  geom_col() +
  coord_flip()+
  scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
diffabund_deseq2 + theme_bw()
ggsave("DESeq2 Barplot.tiff", width = 5, height = 5, device='tiff', dpi=600)

###MaAsLin2 Supplements###
sigtabmaaslin2 <- read.delim2("exp5_pankaj2_output/significant_results.tsv")
diffabund_maaslin2 <- ggplot(sigtabmaaslin2,aes(x=feature,y=coef,fill=coef>0))+
  geom_col() + coord_flip()+
  scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
diffabund_maaslin2 + theme_bw()
ggsave("2023_01_19_exp5_1_diffabund_maaslin2.tiff", width = 5, height = 5, device='tiff', dpi=600)


########SpiecEasi########

#Filter to only OTUs found in 50% of samples and not listed as NA at the Genus level
exp5_pankaj2_prev <- filter_taxa(exp5_pankaj2 , function(x) sum(x >= 1) > (0.50*length(x)), TRUE)
exp5_pankaj2_prev_narm = subset_taxa(exp5_pankaj2_prev, Genus!="NA")

###Pull out the OTU, taxonomy, and metadata tables as data frames###
exp5_pankaj2_narm_otu <- as.data.frame(exp5_pankaj2_prev_narm@otu_table)
exp5_pankaj2_narm_tax <- as.data.frame(exp5_pankaj2_prev_narm@tax_table)
exp5_pankaj2_narm_meta <- as.data.frame(exp5_pankaj2_prev_narm@sam_data)

###Separate out HO and NO###
exp5_pankaj2_narm_otu_ho <- exp5_pankaj2_narm_otu[,which(exp5_pankaj2_narm_meta$Oxygen == "HO")]
exp5_pankaj2_narm_otu_NO <- exp5_pankaj2_narm_otu[,which(exp5_pankaj2_narm_meta$Oxygen == "NO")]

###Prepare HO samples###
exp5_pankaj2_narm_tax<-as.matrix(exp5_pankaj2_narm_tax)
OTU_pankaj2_narm_ho = otu_table(exp5_pankaj2_narm_otu_ho, taxa_are_rows = TRUE)
TAX_pankaj2_narm_ho= tax_table(exp5_pankaj2_narm_tax)
samples_pankaj2_narm_ho= sample_data(exp5_pankaj2_narm_meta)

exp5_pankaj2_narm_ho <- phyloseq(OTU_pankaj2_narm_ho, TAX_pankaj2_narm_ho, samples_pankaj2_narm_ho)

###Prepare NO Samples###
exp5_pankaj2_narm_tax<- as.matrix(exp5_pankaj2_narm_tax)
OTU_pankaj2_narm_NO = otu_table(exp5_pankaj2_narm_otu_NO, taxa_are_rows = TRUE)
TAX_pankaj2_narm_NO= tax_table(exp5_pankaj2_narm_tax)
samples_pankaj2_narm_NO= sample_data(exp5_pankaj2_narm_meta)

exp5_pankaj2_narm_NO <- phyloseq(OTU_pankaj2_narm_NO, TAX_pankaj2_narm_NO, samples_pankaj2_narm_NO)

###Run SPIEC-EASI, set threshold at 0.1 as recommended###                        
se.exp5.pankaj2.ho <- spiec.easi(exp5_pankaj2_narm_ho, method='mb', nlambda=99,
                                lambda.min.ratio=1e-3, pulsar.params = list(thresh = 0.1))

se.exp5.pankaj2.no <- spiec.easi(exp5_pankaj2_narm_NO, method='mb', nlambda=99,
                                lambda.min.ratio=1e-3, pulsar.params = list(thresh = 0.1))

###HO Plot###
spiec.graph.pankaj2.ho=adj2igraph(getRefit(se.exp5.pankaj2.ho), vertex.attr=list(name=taxa_names(exp5_pankaj2_narm_ho)))

narm_spiec_plot_pankaj2_ho <- plot_network(spiec.graph.pankaj2.ho, exp5_pankaj2_narm_ho, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_pankaj2_ho + geom_point(shape = 1,size = 4,colour = "black")
ggsave("SpiecEasi HO.png", width = 15, height = 5, device='png', dpi=600)

narm_spiec_plot_pankaj2_ho_nolegend <- plot_network(spiec.graph.pankaj2.ho, exp5_pankaj2_narm_ho, type='taxa', color="Genus", label = NULL)+
  theme(legend.position = "none")
narm_spiec_plot_pankaj2_ho_nolegend + geom_point(shape = 1,size = 4,colour = "black")
ggsave("SpiecEasi HO no legend.png", width = 10, height = 5, device='png', dpi=600)

#How many OTUs are present
nodes_ho <- gorder(spiec.graph.pankaj2.ho)
#Number of connection between OTUs
edges_ho <- gsize(spiec.graph.pankaj2.ho)
#Number of shortest paths through node, measure of centNOlity
betweenness_ho <- as.list(betweenness(spiec.graph.pankaj2.ho, normalized = TRUE))
#Number of edges per node
degree_ho <- as.data.frame(degree(spiec.graph.pankaj2.ho))

###NO Plot###
spiec.graph.pankaj2.NO=adj2igraph(getRefit(se.exp5.pankaj2.no), vertex.attr=list(name=taxa_names(exp5_pankaj2_narm_NO)))

narm_spiec_plot_pankaj2_NO <- plot_network(spiec.graph.pankaj2.NO, exp5_pankaj2_narm_NO, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_pankaj2_NO + geom_point(shape = 1,size = 4,colour = "black")
ggsave("SpiecEasi NO.png", width = 15, height = 5, device='png', dpi=600)

narm_spiec_plot_pankaj2_NO_nolegend <- plot_network(spiec.graph.pankaj2.NO, exp5_pankaj2_narm_NO, type='taxa', color="Genus", label = NULL)+
  theme(legend.position = "none")
narm_spiec_plot_pankaj2_NO_nolegend + geom_point(shape = 1,size = 4,colour = "black")
ggsave("SpiecEasi NO no legend.png", width = 10, height = 5, device='png', dpi=600)

nodes_NO <- gorder(spiec.graph.pankaj2.NO)
edges_NO <- gsize(spiec.graph.pankaj2.NO)
betweenness_NO <- as.list(betweenness(spiec.graph.pankaj2.NO, normalized = TRUE))
degree_NO <- as.data.frame(degree(spiec.graph.pankaj2.NO))

###Degree plot###
degree_df <- data.frame(degree_ho,degree_NO)
colnames(degree_df) <- c("HO","NO")

count_ho <- as.data.frame(table(degree_df$HO))
colnames(count_ho) <- c("Degree","HO")

count_NO <- as.data.frame(table(degree_df$NO))
colnames(count_NO) <- c("Degree","NO")

degree_freq <- rbind.fill(count_ho[c("Degree", "HO")], count_NO[c("Degree", "NO")])

degree_freq <- degree_freq %>%
  group_by(Degree) %>%
  dplyr::summarise(across(c(HO,NO), ~sum(., na.rm = TRUE)))

degree_long <- degree_freq %>%
  pivot_longer(c("HO","NO"), names_to = "Oxygen", values_to = "Frequency")

#Plot how many nodes have a certain degree for both groups
ggplot(degree_long,aes(x=Degree,y=Frequency,group=Oxygen))+
  geom_line(aes(color=Oxygen)) +
  geom_point(aes(color=Oxygen)) +
  theme_bw()
ggsave("Degree Plot.png", width = 5, height = 5, device='png', dpi=600)

###Edges Plot###
edges_df <- data.frame(edges_ho,edges_NO)
colnames(edges_df) <- c("HO","NO")
edges_long <- edges_df %>%
  pivot_longer(c("HO","NO"), names_to = "Oxygen", values_to = "Edges")

#Plot total edges in each group
ggplot(edges_long,aes(x= Oxygen, y = Edges))+
  geom_bar(stat="identity", aes(fill=Oxygen))+
  theme_bw()
ggsave("Edge Plot.png", width = 5, height = 5, device='png', dpi=600) 


########Random Forest########

exp5_combined_gen_prev_noNA = subset_taxa(exp5_pankaj2_gen_prev, Genus!="NA")
exp5_combined_gen_prev_rel <- transform_sample_counts(exp5_combined_gen_prev_noNA, function(x) x / sum(x) )

###Prep OTU table###
predictors <- as.data.frame(exp5_combined_gen_prev_rel@otu_table)
taxlabels <- as.data.frame(exp5_combined_gen_prev_rel@tax_table)
predictors <- t(predictors)
colnames(predictors) <- taxlabels$Genus
#Get relevant metadata
response <- as.factor(sample_data(exp5_combined_gen_prev_rel)$Oxygen)
alsopredictors <- sample_data(exp5_combined_gen_prev_rel)[,4:7]
#Surprise tool that will help us later
tax.labels <- as.data.frame(exp5_combined_gen_prev_rel@tax_table)

###Make dataframe combining OTUs and Metadata###
rf.data <- data.frame(response, predictors, alsopredictors)

###Set random seed and partition data into training and testing###
set.seed(1312)
pankaj2_idx = createDataPartition(rf.data$response, p = 0.75, list = FALSE)
pankaj2_train = rf.data[pankaj2_idx, ]
pankaj2_test = rf.data[-pankaj2_idx, ]

###Set cross-validation, in this case Leave One Out Cross-Validation###
pankaj2_cv <- trainControl(method='LOOCV', classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = TRUE)

#Weight model so both groups are weighted equally
model_weights <- ifelse(pankaj2_train$response == "HO",
                        (1/table(pankaj2_train$response)[1]) * 0.5,
                        (1/table(pankaj2_train$response)[2]) * 0.5
)

###Set paNOmeter tuning grid, testing a few paNOmeter combinations###
tgrid <- expand.grid(
  .mtry = c(2,16,32,63),
  .splitrule = c("gini","extratrees"),
  .min.node.size = c(1,3,5,10,18)
)

###Run random forest###
set.seed(1312)
rf_fit <- train(as.factor(response) ~ ., 
                data = pankaj2_train, 
                method = "ranger",
                tuneGrid = tgrid,
                num.trees = 1000,
                trControl = pankaj2_cv,
                weights = model_weights,
                metric = "ROC",
                verbose = FALSE,
                verbosity = 0,
                importance = "permutation"
)

rf_fit

###New grid, using the best hyperpaNOmeter combination###
best_grid <- expand.grid(
  .mtry = 63,
  .splitrule = "gini",
  .min.node.size = 1
)

###Final random forest model###
set.seed(1312)
rf_fit2 <- train(as.factor(response) ~ ., 
                 data = pankaj2_train, 
                 method = "ranger",
                 tuneGrid = best_grid,
                 num.trees = 1000,
                 trControl = pankaj2_cv,
                 weights = model_weights,
                 metric = "Sens",
                 verbose = FALSE,
                 verbosity = 0,
                 importance = "permutation"
)

rf_fit2

###Plot ROC###
res_rf <- evalm(rf_fit2)
ggsave("Random Forest ROC.png", width = 5, height = 5, device='png', dpi=600)

###Check against test data###
testclass <- predict(rf_fit2, newdata = pankaj2_test)
cfMatrix <- confusionMatrix(data = testclass, pankaj2_test$response)

###Plot feature importance###
features <- varImp(rf_fit2)
features <- features$importance
features <- rownames_to_column(features)
features <- features[order(-features$Overall),]
features_top50 <- as.data.frame(features[1:50,])
features_top25 <- as.data.frame(features[1:25,])

ggplot2::ggplot(features_top50, aes(x=rowname, y=Overall)) +
  geom_point( color="#FB0207", size=4, alpha=1.0)+
  geom_segment( aes(x=rowname, xend=rowname, y=0, yend=Overall), 
                color='#EF94A2') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light() +
  coord_flip() 

ggsave("Random Forest Top Features.png", width = 5, height = 10, device='png', dpi=600)

######Spearman Correlation of Important OTUs with qPCR Data#####
exp5_pankaj2_meta_df <- as.data.frame(exp5_pankaj2_meta)
labeled_otus <- exp5_pankaj2_otu
colnames(labeled_otus) <- exp5_pankaj2_tax$Genus
mv.combined.df <- data.frame(exp5_pankaj2_meta_df,labeled_otus)
cor_df <- mv.combined.df[,4:228]
cor_df <- cor_df[,c(4,3,2,1,5,6,13,20,22,50,62,119)]


tiff("Spearman Correlation Significant3.tif", width=5*dpi, height=5*dpi, res=dpi)
corr<-cor(cor_df[,5:ncol(cor_df)], cor_df[,1:4],method = "spearman")
cor_test_mat <- corr.test(cor_df[,5:ncol(cor_df)], cor_df[,1:4], method = "spearman",adjust = "fdr")$p
ggcorrplot(corr, lab = TRUE, lab_size = 3)
dev.off()

