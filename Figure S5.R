library(tidyverse)
library(vegan)
library(phyloseq)
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
library(psych)
library(ggcorrplot)

########################## Figure S5B ##########################

########### Renyi Diversity Profile ###########

###Get Renyi diversity profile using Hill numbers 1-4###
set.seed(420)
renyihill <- renyi(pankaj2.sort.flip, scales = c(1, 2, 3, 4),hill = TRUE)
renyihill2 <- renyihill %>%
  mutate(
    Oxygen = md.pankaj2.sort$Oxygen
  )
renyihill2 <- rownames_to_column(renyihill2, var = "Name")
renyihill2_long <- renyihill2 %>% pivot_longer(cols=c(2:5),
                                               names_to='Hill',
                                               values_to='Value')

###Plot line plot###
ggplot(renyihill2_long,aes(x=Hill,y=Value,group=Name))+
  geom_line(aes(color=Oxygen)) +
  geom_point(aes(color=Oxygen), size = 3) +
  scale_color_manual(values = deseq2_colors) +
  theme_bw()
#Save

########### Weighted UniFrac ###########

###Load Phylogenetic tree###
tree <- ape::read.nexus("AMP_All_Tree.nxs")

###Fix labels in the phylogenetic tree and subset taxonomy to what is found in the tree###
labels <- as.data.frame(tree[["tip.label"]])
labels[c('Taxa', 'Name')] <- str_split_fixed(labels$`tree[["tip.label"]]`, ";", 2)
newlabel <- gsub("(^')|('$)", "", labels$Name)
tree[["tip.label"]] <- newlabel

pankaj2.clean.new <- pankaj2.clean[row.names(pankaj2.clean) %in% newlabel,]
pankaj2.taxa.new <- pankaj2.taxa[row.names(pankaj2.taxa) %in% newlabel,]

###Create new phyloseq object with tree###
OTU.with.tree = otu_table(pankaj2.clean.new, taxa_are_rows = TRUE)
TAX.with.tree = tax_table(pankaj2.taxa.new)
samples.with.tree = sample_data(md.pankaj2.sort2)
phylo <- phy_tree(tree)

exp5_pankaj2_withtree <- phyloseq(OTU, TAX, samples, phylo)
exp5_pankaj2_gen_withtree <- tax_glom(exp5_pankaj2_withtree, taxrank = 'Genus')

###Plot phylogenetic tree###
plot_tree(exp5_pankaj2_gen_withtree, nodelabf=nodeplotblank, label.tips="Genus", ladderize="left", color = "Phylum")
#Save

###Get UniFrac Distance###

exp5_pankaj2_rel_withtree <- transform_sample_counts(exp5_pankaj2_gen_withtree, function(x) x / sum(x) )

uni <- phyloseq::UniFrac(exp5_pankaj2_rel_withtree, weighted = TRUE)
ord = ordinate(exp5_pankaj2_rel_withtree, "PCoA", "unifrac", weighted=TRUE)

exp5_pankaj2_meta_withtree <- as.data.frame(exp5_pankaj2_rel_withtree@sam_data)
exp5_pankaj2_phylo <- phy_tree(exp5_pankaj2_rel_withtree)

factor_exp5_unifrac <- as.factor(exp5_pankaj2_meta_withtree$Oxygen)
type_exp5_unifrac <- as.numeric(factor_exp5_unifrac)
pca_colors <- c("#FB0207","#EF94A2")

#Plot
dpi=600
tiff("Figure S5B UniFrac.tif", width=5*dpi, height=5*dpi, res=dpi)
plot(c(-0.2, 0.6), c(-0.2, 0.15), font = 2, font.lab = 2, xlab="PC1", ylab="PC2", type="n")
points(ord$vectors[,1:2], pch = 21, cex = 1.3, bg = pca_colors[type_exp5_unifrac], lwd = 1)
ordiellipse(ord$vectors[,1:2], factor_exp5_unifrac, kind = "se", conf = 0.95,col = pca_colors)
ordispider(ord$vectors[,1:2], factor_exp5_unifrac, label = TRUE)
dev.off()

###PERMANOVA and PERMDISP###
permanova_pcoa_df <- data.frame(exp5_pankaj2_meta)
set.seed(1312)
permanova_pankaj2_pcoa <- vegan::adonis2(uni ~ Oxygen, data = permanova_pcoa_df, permutations = 10000)

disp_pankaj2_pcoa <- betadisper(uni, permanova_pcoa_df$Oxygen)
set.seed(1312)
permdisp_pankaj2_pcoa <- permutest(disp_pankaj2_pcoa, permutations = 10000)

print(permanova_pankaj2_pcoa)
print(permdisp_pankaj2_pcoa)


########################## Figure S5C ##########################

###Get the significant OTUs from the DESeq2 data in Figure 1###
diffabund_deseq2 <- ggplot(sigtab_otu,aes(x=Genus,y=log2FoldChange,fill=log2FoldChange>0))+
  geom_col() +
  coord_flip()+
  scale_fill_manual(values=deseq2_colors, labels=c("negative","positive"))
diffabund_deseq2 + theme_bw()
#Save

########################## Figure S5D ##########################

###Remove OTUs unassigned at Genus level###
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

###Set hyperparameter tuning grid, testing a few parameter combinations###
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

###New grid, using the best hyperparameter combination###
best_grid <- expand.grid(
  .mtry = 63,
  .splitrule = "gini",
  .min.node.size = 5
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

########### ROC ###########
res_rf <- evalm(rf_fit2)
#Save

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

########### Feature Importance Plot ###########
ggplot2::ggplot(features_top50, aes(x=rowname, y=Overall)) +
  geom_point( color="#FB0207", size=4, alpha=1.0)+
  geom_segment( aes(x=rowname, xend=rowname, y=0, yend=Overall), 
                color='#EF94A2') +
  xlab('Variable')+
  ylab('Overall Importance')+
  theme_light() +
  coord_flip() 
#Save

########################## Figure S5E ##########################

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


######## HO Network Plot ########

spiec.graph.pankaj2.ho=adj2igraph(getRefit(se.exp5.pankaj2.ho), vertex.attr=list(name=taxa_names(exp5_pankaj2_narm_ho)))

narm_spiec_plot_pankaj2_ho <- plot_network(spiec.graph.pankaj2.ho, exp5_pankaj2_narm_ho, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_pankaj2_ho + geom_point(shape = 1,size = 4,colour = "black")
#Save

narm_spiec_plot_pankaj2_ho_nolegend <- plot_network(spiec.graph.pankaj2.ho, exp5_pankaj2_narm_ho, type='taxa', color="Genus", label = NULL)+
  theme(legend.position = "none")
narm_spiec_plot_pankaj2_ho_nolegend + geom_point(shape = 1,size = 4,colour = "black")
#Save


######## NO Network Plot ########

spiec.graph.pankaj2.NO=adj2igraph(getRefit(se.exp5.pankaj2.no), vertex.attr=list(name=taxa_names(exp5_pankaj2_narm_NO)))

narm_spiec_plot_pankaj2_NO <- plot_network(spiec.graph.pankaj2.NO, exp5_pankaj2_narm_NO, type='taxa', color="Genus", label = NULL)
narm_spiec_plot_pankaj2_NO + geom_point(shape = 1,size = 4,colour = "black")
#Save

narm_spiec_plot_pankaj2_NO_nolegend <- plot_network(spiec.graph.pankaj2.NO, exp5_pankaj2_narm_NO, type='taxa', color="Genus", label = NULL)+
  theme(legend.position = "none")
narm_spiec_plot_pankaj2_NO_nolegend + geom_point(shape = 1,size = 4,colour = "black")
#Save


########################## Figure S5F ##########################

exp5_pankaj2_meta_df <- as.data.frame(exp5_pankaj2_meta)
labeled_otus <- exp5_pankaj2_otu
colnames(labeled_otus) <- exp5_pankaj2_tax$Genus
mv.combined.df <- data.frame(exp5_pankaj2_meta_df[,-c(4)],labeled_otus)
cor_df <- mv.combined.df[,c(4:7,12,39,93,152,202,207,217)]


tiff("Figure S5F.tif", width=5*dpi, height=5*dpi, res=dpi)
corr<-cor(cor_df[,5:ncol(cor_df)], cor_df[,1:4],method = "spearman")
cor_test_mat <- corr.test(cor_df[,5:ncol(cor_df)], cor_df[,1:4], method = "spearman",adjust = "fdr")$p
ggcorrplot(corr, lab = TRUE, lab_size = 3)
dev.off()
