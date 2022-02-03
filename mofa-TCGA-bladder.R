#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("MOFA2")
#BiocManager::install("MOFAdata")

library(mixOmics)

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)



clean_all_zeros<-function(df){
  df<-as.data.frame(apply(df, 1, function(x) as.numeric(x)))
  ind <- apply(df, 1, function(x) sum(x, na.rm = TRUE)==0) 
  #remove genes with zero variance
  df<-df[!ind,]
  return(df)
}


# depends on SPLS.R to process X1,X2
X1_mat<-clean_all_zeros(X1_raw[,-1])
X2_mat<-clean_all_zeros(X2_raw[,-1])

df<-X1_raw[,-1]


data<-list()
data$mRNA<-as.matrix(t(X1_t))
data$protein<-as.matrix(t(X2_t))
  
MOFAobject <- create_mofa(data)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
data_opts
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 4
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

outfile = file.path(getwd(),"model.hdf5")
MOFAobject <- run_mofa(MOFAobject, outfile)

#MOFAobject<-MOFAobject.trained

#metadata<- data.frame(subtype=breast.TCGA$data.train$subtype)


#metadata<-cbind(MOFAobject@samples_metadata, metadata)

#samples_metadata(MOFAobject)<-metadata

plot_factor_cor(MOFAobject)
saveRDS(MOFAobject, 'bladder_cancer/mofa_object.RDS')
plot_variance_explained(MOFAobject, max_r2=14)


plot_variance_explained(MOFAobject, plot_total = T)[[2]]


plot_factor(MOFAobject.trained, 
            factors = 1, 
            color_by = "Factor1"
)

dev.off()
plot_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject,
                 view = 'protein',
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


plot_top_weights(MOFAobject,
                 view = 'mRNA',
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
output<-'bladder_cancer/plots/mofa/'
ggsave(paste0(output,'top_weights.png'), width = 4, height = 4, dpi=500 )
ggsave(paste0(output,'GSEA_factor_',factor_to_plot,'.png'), width = 9, height=4, dpi=100)

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

factor_to_plot=1
plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = factor_to_plot,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


ggsave(paste0(output,'heatmap',factor_to_plot,'.png'), width = 9, height=4, dpi=100)

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
data("reactomeGS", package = "MOFAdata")

colors <- c("#999999", "#E69F00", "#56B4E9")
colors <- colors[as.numeric(MOFAobject@samples_metadata$subtype)]

MOFAobject@expectations$W


scatterplot3d(MOFAobject@expectations$Z$group1[,c(1,2,3)], pch=16, color=colors)
#install.packages('GGally')

plot_factors(MOFAobject, 
             factors = c(1:3), 
             color_by = 'subtype'
)

#### Enrichment analysis 
#utils::data()

head(colnames(reactomeGS))
require('biomaRt')
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

## Bimap interface:
## select() interface:
## Objects in this package can be accessed using the select() interface
## from the AnnotationDbi package. See ?select for details.
## Bimap interface:
x <- org.Hs.egGENENAME
# Get the gene names that are mapped to an entrez gene identifier
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
if(length(xx) > 0) {
  # Get the GENE NAME for the first five genes
  xx[1:5]
  # Get the first one
  xx[[1]]
}

mart <- useMart('ENSEMBL_MART_ENSEMBL')
mart <- useDataset('hsapiens_gene_ensembl', MOFAobject$m)

annotLookup <- getBM(
  mart = mart,
  attributes = c(
    'hgnc_symbol',
    'ensembl_gene_id',
    'gene_biotype'),
  uniqueRows = TRUE)



# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)

# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
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

res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)

library("AnnotationDbi")
install.packages('org.Hs.eg.db')
library("org.Hs.eg.db")
BiocManager::install("org.Hs.eg.db")


ensid = mapIds(org.Hs.eg.db,
                  keys=rownames(MOFAobject@data$mRNA$group1), 
                  column="ENSEMBL",
                  keytype="SYMBOL",
                  multiVals="first")



select(edb, keys=keys, columns=rownames(MOFAobject@data$mRNA$group1),
       keytype="GENEID")

