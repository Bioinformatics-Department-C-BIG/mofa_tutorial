reticulate::use_python(python = "C:/Users/athienitie/Anaconda3/python.exe")


library(mixOmics)
data(breast.TCGA)
#BiocManager::install("biomaRt")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)


breast_data<-list()
breast_data$miRNA<-t(as.matrix(breast.TCGA$data.train$mirna))
breast_data$mRNA<-t(as.matrix(breast.TCGA$data.train$mrna))
breast_data$protein<-t(as.matrix(breast.TCGA$data.train$protein))



  
MOFAobject <- create_mofa(breast_data)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 13
model_opts



train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

train_opts


MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)


MOFAobject <- run_mofa(MOFAobject, outfile="MOFA2_CLL.hdf5")


metadata<- data.frame(subtype=breast.TCGA$data.train$subtype)


metadata<-cbind(MOFAobject@samples_metadata, metadata)

samples_metadata(MOFAobject)<-metadata

plot_factor_cor(MOFAobject)

plot_variance_explained(MOFAobject, max_r2=14)


plot_variance_explained(MOFAobject, plot_total = T)[[2]]


plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)


plot_weights(MOFAobject,
             view = "mRNA",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject,
                 view = 'protein',
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)



#####Prediction of clinical subgroups 

suppressPackageStartupMessages(library(randomForest))


# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

# Train the model for IGHV
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
model.ighv <- randomForest(IGHV ~ ., data=df[!is.na(df$IGHV),], ntree=10)
df$IGHV <- NULL

# Do predictions
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)



  
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "ER-alpha",
            add_violin = TRUE,
            dodge = TRUE
)



plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "ER-alpha"
) + labs(y="RNA expression")


plot_data_heatmap(MOFAobject, 
                  view = "miRNA",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)



p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "subtype",
                  shape_by = "subtype",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

###### TODO: overlay the subtypes in 3d 
#install.packages("scatterplot3d") # Install
library("scatterplot3d") # load


colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(MOFAobject@samples_metadata$subtype)]

MOFAobject@expectations$W


scatterplot3d(MOFAobject@expectations$Z$group1[,c(1,2,3)], pch=16, color=colors)
#install.packages('GGally')

library()
plot_factors(MOFAobject, 
             factors = c(1:3), 
             color_by = 'subtype'
)

#### Enrichment analysis 
utils::data()

head(colnames(reactomeGS))


# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reac, 
                               view = "mRNA",
                               sign = "negative"
)

#install.packages("plot3D")
###Convert gene names #
library("plot3D")
library("tidyr")
library (biomaRt)

#separate(MOFAobject@expectations$Z$group1[,c(1,2,3)])
#scatter3D(, colvar = NULL, col = "blue",
#          pch = 19, cex = 0.5)
install.packages(b)
library (biomaRt)

g_All <- getGene(id = rownames(breast_data$mRNA) , type='hgnc_symbol' ,mart=ensembl )


