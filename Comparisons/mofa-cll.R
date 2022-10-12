#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install('EnsDb.Hsapiens.v79')
library(EnsDb.Hsapiens.v79)

map_to_gene<-function(ens_ids){
  ensembldb::select(EnsDb.Hsapiens.v79, keys= ens_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
}

#BiocManager::install("MOFA2")
#devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
#browseVignettes("MOFA2")
#BiocManager::install("MOFAdata")
#BiocManager::install('MultiAssayExperiment')

outdir<-'Comparisons/plots/mofa/'

library(MultiAssayExperiment)
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)


utils::data("CLL_data")  
data("CLL_data")



lapply(CLL_data,dim)

CLL_data$Drugs

CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")

#RENAME TO ENSids 



MOFAobject <- create_mofa(CLL_data)
plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10

model_opts$likelihoods


#### Training 
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

train_opts


##### Prepare MOFA 
MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts
)

##### slots 
#MOFAobject <- MOFA2::run_mofa(MOFAobject)

MOFAobject <- readRDS(url("http://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/MOFA2_CLL.rds"))

slotNames(MOFAobject)
names(MOFAobject@data)

dim(MOFAobject@data$Drugs$group1)

names(MOFAobject@expectations)

dim(MOFAobject@expectations$Z$group1)

# mrna weight matrix 
dim(MOFAobject@expectations$W$mRNA)

MOFAobject@features_metadata$feature
## change names 

samples_metadata(MOFAobject) <- CLL_metadata

MOFAobject_gs<-MOFAobject
MOFA2::features_names(MOFAobject)

ens_ids<- features_names(MOFAobject)$mRNA
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ens_ids, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
features_names(MOFAobject_gs)$mRNA<-geneIDs1$SYMBOL
plot_factor_cor(MOFAobject)




# Result 1: Look at all the variances

plot_variance_explained(MOFAobject, max_r2=15, factors = 1:9)
ggsave(paste0(outdir, 'variance_explained','.png'), width = 4, height=4, dpi=100)



TOP<-MOFA2::get_weights(MOFAobject,views = 'mRNA', factors=1, 
                   as.data.frame =TRUE)


TOP[1:100,]
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


#### Association 
### Test the association between MOFA factors and gender, survival outcome 
### and age 

correlate_factors_with_covariates(MOFAobject, 
                                  covariates = c("Gender","died","age"), 
                                  plot="log_pval"
)


# todo get top variables!! 


#### Factor analysis 
### like principal component analysis 

  plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)



### Plot feature weights 
### Feature weights for mutations 
fps<-seq(1:5)

vps<-c("Mutations", "Drugs", 'Methylation', "mRNA")

v_set=c()
for (i in 1:length(vps)){
  for (ii in 1:length(fps)){

plot_top_weights(MOFAobject_gs,
             view = vps[i],
             factor = fps[ii],
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
ggsave(paste0(outdir, 'top_weights_', fps[ii],'_',vps[i],'.png'), width = , height=4, dpi=100)

plot_weights(MOFAobject_gs, 
             view = vps[i], 
             factor = fps[ii], 
             nfeatures = 10
)
ggsave(paste0(outdir, 'all_weights_', fps[ii],'_',vps[i],'.png'), width = 4, height=4, dpi=100)

  
  
  ###### Heatmaps 
nfs=20

jpeg(paste0(outdir, 'heatmap_',fps[ii],'_','mRNA_', 'nfs_', nfs, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
fps[ii]=1
# Plot heatmaps for each factor only for mRNA 
p<-plot_data_heatmap(MOFAobject_gs, 
                    view = 'mRNA', 
                    factor =  fps[ii],  
                    features = nfs,
                    denoise = TRUE,
                    cluster_rows = FALSE, cluster_cols = FALSE,
                    show_rownames = TRUE, show_colnames = FALSE,
                    scale = "row",
                    fontsize_number = 5
                    
                    
)
dev.off()
  
  jpeg(paste0(outdir, 'heatmap_',fps[ii],'_','Drugs_', 'nfs_', nfs, '.jpeg'), height=20*nfs, width=20*nfs)
  
  p<-plot_data_heatmap(MOFAobject_gs, 
                       view = 'Drugs', 
                       factor =  fps[ii],  
                       features = nfs,
                       denoise = TRUE,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       show_rownames = TRUE, show_colnames = FALSE,
                       scale = "row"
  )
  dev.off()
  
  # top weights
  # concat all 

 
  
  }
  
 
}


#one eheatmap
for (ii in 1:5){
jpeg(paste0(outdir, 'heatmap_',ii,'_','Methylation_', 'nfs_', nfs, '.jpeg'), height=20*nfs, width=20*nfs)
p<-plot_data_heatmap(MOFAobject_gs, 
                     view = 'Methylation', 
                     factor =  ii,  
                     features = nfs,
                     denoise = TRUE,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     show_rownames = TRUE, show_colnames = FALSE,
                     scale = "row"
)
dev.off()
}

v_set=c()
for (i in 1:15){
  weights<-MOFA2::get_weights(MOFAobject_gs,views = 'mRNA', factors = i)
  top_weights<-rownames(weights$mRNA)[order(abs(weights$mRNA), decreasing = TRUE) ][1:15]
  #top_genes_to_compare<-map_to_gene(top_weights)$SYMBOL
 top_genes_to_compare<-top_weights
  write.csv(top_genes_to_compare, paste0(outdir, 'top_genes_', i, '_', '.csv'))
  print(top_genes_to_compare)
  # concat all 
  v_set <- c(v_set, top_genes_to_compare)
}

write.csv(v_set, paste0(outdir, 'all_genes_top_', total_comps, '.csv'))

###Plot gene weights for MRNA expression #


MOFA2::features_names(MOFAobject)$mRNA

map_names<-rownames(MOFAobject@data$mRNA$group1)



#   comparison no1: top 30 genes 
fp<-1




jpeg(paste0(outdir, 'mofa_factor_space', 1,'.jpeg'))

#VECTOR of 2 groups 
combined_groups<-factor(paste(MOFAobject@samples_metadata$IGHV,MOFAobject@samples_metadata$trisomy12, sep="_"))
notnas<-combined_groups %in% c('0_0', '0_1', '1_0', '1_1')
combined_groups[!notnas]<-NA
combined_groups<-factor(combined_groups)
levels(combined_groups)<-c('IGHV-,no tr12', 'IGHV-,tr12', 'IGHV+,no tr12', 'IGHV+,tr12')
### experiment 1 show ighv and trisomy 12

fps_sel=c(1,3)

jpeg(paste0(outdir, 'mofa_factor_space_new', paste0(fps_sel,collapse='_'),'.jpeg'), res=150, width = 6, height=5, units='in')



p<-plot_factors(MOFAobject,
                  factors = fps_sel, 
                  dot_size = 2,
                  color_name = 'IGHV',
                  color_by=combined_groups, 
                  show_missing = T)
p+theme(text=element_text(size=14))

dev.off()
utils::data(reactomeGS)



##### Inspection of combinations of Factors

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_by = "IGHV",
                  shape_by = "trisomy12",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)



p
dev.off()

##### Predicton of cilincal subgroups 
install.packages('randomForest')
suppressPackageStartupMessages(library(randomForest))


df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

df <- as.data.frame(get_factors(MOFAobject, factors=c(1,2))[[1]])

# Train the model for IGHV
df$IGHV <- as.factor(MOFAobject@samples_metadata$IGHV)
model.ighv <- randomForest(IGHV ~ ., data=df[!is.na(df$IGHV),], ntree=10)
df$IGHV <- NULL

# Do predictions
MOFAobject@samples_metadata$IGHV.pred <- stats::predict(model.ighv, df)



# Train the model for Trisomy12
df$trisomy12 <- as.factor(MOFAobject@samples_metadata$trisomy12)
model.trisomy12 <- randomForest(trisomy12 ~ ., data=df[!is.na(df$trisomy12),], ntree=10)
df$trisomy12 <- NULL

MOFAobject@samples_metadata$trisomy12.pred <- stats::predict(model.trisomy12, df)


MOFAobject@samples_metadata$IGHV.pred_logical <- c("True","Predicted")[as.numeric(is.na(MOFAobject@samples_metadata$IGHV))+1]

p <- plot_factors(MOFAobject, 
                  factors = c(1,3), 
                  color_name = "IGHV.pred",
                  shape_name = "IGHV.pred_logical",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)


###### GSEA 
library(MOFAdata)

utils::data(reactomeGS)

head((reactomeGS))

head(colnames(reactomeGS))

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

for (ii in 1:10){
  sign<-res.positive$pval.adj[,ii][which(res.positive$pval.adj[,ii]<0.05)]
  #sign<-res.negative$pval.adj[,ii][which(res.negative$pval.adj[,ii]<0.005)]
  
  print(length(sign))
}
MOFAobject@data$Mutations
names(res.positive)[which(res.positive$pval.adj<0.05)]

View(res.positive$pval.adj)


theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

names(res.positive)

# ENS ids seems to be required for the enrichment 
plot_enrichment_heatmap(res.positive)
ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)
dev.off()

plot_enrichment_heatmap(res.negative)
ggsave(paste0(outdir,'Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=120)


factor_to_plot=5
nfs=7
plot_enrichment(res.positive, factor = factor_to_plot, 
                max.pathways = nfs)+
  theme(text=element_text(size=12))
ggsave(paste0(outdir,'GSEA_factor_',factor_to_plot, '_', nfs, '.jpeg'), width = 5, height=4, dpi=200)


plot_enrichment(res.negative, factor = factor_to_plot, max.pathways = nfs)
ggsave(paste0(outdir,'GSEA_factor_neg_',factor_to_plot, '_', nfs,'.jpeg'), width = 5, height=4, dpi=200)


##### ASSOCIATE FACTORS WITH SURVIVAL 
library(survival)
library(survminer)


SurvObject <- Surv(MOFAobject@samples_metadata$TTT, MOFAobject@samples_metadata$treatedAfter)
Z <- get_factors(MOFAobject)[[1]]
fit <- coxph(SurvObject ~ Z) 
fit


##### Plot hazard ratios 
s <- summary(fit)
coef <- s[["coefficients"]]

df <- data.frame(
  factor = factor(rownames(coef), levels = rev(rownames(coef))),
  p      = coef[,"Pr(>|z|)"], 
  coef   = coef[,"exp(coef)"], 
  lower  = s[["conf.int"]][,"lower .95"], 
  higher = s[["conf.int"]][,"upper .95"]
)

jpeg(paste0(outdir, 'hazard_ratio_vs_factors_', 1,'_','mRNA','.jpeg'), res=120)

options(digits=2)
ggplot(df[1:10,], aes(x=factor, y=coef, ymin=lower, ymax=higher, label= format(p, scientific = TRUE))) +
 geom_pointrange( col='#619CFF') +
  geom_text(hjust=-1.7)+
  coord_flip() +
  scale_x_discrete() + 
  labs(y="Hazard Ratio", x="") + 
  geom_hline(aes(yintercept=1), linetype="dotted")+
  theme(text=element_text( size=17),
        plot.margin=margin(1,2,1,2,unit='cm') )
dev.off()


#### MAP TO GENE SYMBOLS FOR PLOTTING 
ensemblsIDS<-rownames(MOFAobject@data$mRNA$group1)
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")

symbols

#ANOTHER SOLUTION 
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
symbols2<-getBM(attributes='hgnc_symbol', 
      filters = 'ensembl_gene_id', 
      values = ensemblsIDS, 
      mart = ensembl)

