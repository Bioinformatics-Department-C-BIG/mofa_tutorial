
#install.packages('psych')


##### 



plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained','.png'), width = 4, height=4, dpi=100)


#MOFAobject@samples_metadata$Grade<-as.factor(MOFAobject@samples_metadata$Grade)
MOFAobject@samples_metadata$Grade
colnames(MOFAobject@samples_metadata)


stats<-apply(MOFAobject@samples_metadata, 2,table )

#### Covariance of factors with metadata 
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c('INFODT.x','NUPSOURC',"NP1ANXS", "NP1DPRS", 'NP1HALL', 'NP1APAT', 'NP1DDS', 'NP1RTOT'), 
                                  plot = "log_pval",
                                  
)

correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("NP1ANXS", "NP1DPRS", 'NP1HALL', 'NP1APAT', 'NP1DDS', 'NP1RTOT'), 
                                  plot = "log_pval",
                                  
)
dev.off()

correlate_factors_with_covariates(MOFAobject,
                                  covariates = c('INFODT.x','SCAU1','SCAU2',"SCAU3", "SCAU4", "SCAU5", "SCAU6",  "SCAU7", "SCAU8", 'NP1HALL', 'NP1APAT', 'NP1DDS', 'NP1RTOT'), 
                                  plot = "log_pval",
                                  
)
ggsave(paste0(outdir, 'factors_covariates','.png'), width = 4, height=4, dpi=100)


view='miRNA'; factor=2
all_weights<-MOFA2::get_weights(MOFAobject,views = view, factors=factor, 
                        as.data.frame =TRUE)

### get the top highly weighted variables - absolute value
top<-all_weights[order(abs(all_weights$value), decreasing = TRUE),][1:10,]
write.table(top,paste0(outdir, 'top_weights_vals',factor,'_', view,'.txt'), sep = '\t')

   plot_variance_explained(MOFAobject, plot_total = T)[[2]]
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


#install.packages('psych')

plot_factors(MOFAobject, 
             factors = c(1,2), 
             dot_size = 2.5
)

plot_weights(MOFAobject,
             view = "miRNA",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)





### Age, gender , stage does not discriminate factors
#Conclusion here:# factor 1 correlates with grade

p<-plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "EORTC.risk",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)
p
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Subtype",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Subtype",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)

# Factor 2 associates with proteomic Subtype 

plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Subtype",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)

##### plot weights 
fps<-seq(1:5)

vps<-c("proteomics" ,"miRNA")

v_set=c()
v_set=c()
fps<-seq(1:3)

view='miRNA'
factor=3
plot_top_weights(MOFAobject,
                 view = view,
                 factor = factor,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
dev.off()
ggsave(paste0(outdir, 'top_weights_',factor, view,'_','.png'), width =3 , height=4, dpi=100)



for (i in 1:length(vps)){
  for (ii in 1:length(fps)){
    
    plot_top_weights(MOFAobject,
                     view = vps[i],
                     factor = fps[ii],
                     nfeatures = 10,     # Top number of features to highlight
                     scale = T           # Scale weights from -1 to 1
    )
    ggsave(paste0(outdir, 'top_weights_', fps[ii],'_',vps[i],'.png'), width = , height=4, dpi=100)
    
    plot_weights(MOFAobject, 
                 view = vps[i], 
                 factor = fps[ii], 
                 nfeatures = 10
    )
    ggsave(paste0(outdir, 'all_weights_', fps[ii],'_',vps[i],'.png'), width = 4, height=4, dpi=100)
    
    
    
    ###### Heatmaps 
    nfs=20
    
    jpeg(paste0(outdir, 'heatmap_',fps[ii],'_','miRNA_', 'nfs_', nfs, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
    fps[ii]=1
    # Plot heatmaps for each factor only for miRNA 
    p<-plot_data_heatmap(MOFAobject, 
                         view = vps[i], 
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
    
    p<-plot_data_heatmap(MOFAobject, 
                         view = vps[i], 
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








plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 1,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "miRNA",
                  factor = 1,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "miRNA",
                  factor = 3,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 3,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


### So multi omics factors are more related to Stage than to subtype!!
##How?? plot on the factors and color by stage

p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "EORTC.risk",
                  shape_by = "Grade",
                  dot_size = 2.5,
                  show_missing = T
)

#p <- p + 
#  geom_hline(yintercept=-1, linetype="dashed") +
#  geom_vline(xintercept=(-0.5), linetype="dashed")
print(p)
ggsave(p)
ggsave(paste0(outdir,'factor_plot','.png'), width = 4, height=4, dpi=120)





#### Now make predictions
### 



### 
## GENE SET ENRICHMENT! 
## AND reactome gsa enrichment!!


###### GSEA 
library(MOFAdata)

utils::data(reactomeGS)

head((reactomeGS))

features_names(MOFAobject)$miRNA
# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "miRNA",
                               sign = "positive"
)



# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "miRNA",
                               sign = "negative"
)

for (ii in 1:4){
  sign<-res.positive$pval.adj[,ii][which(res.positive$pval.adj[,ii]<0.05)]
  #sign<-res.negative$pval.adj[,ii][which(res.negative$pval.adj[,ii]<0.005)]
  
  print(length(sign))
}
names(res.positive)[which(res.positive$pval.adj<0.05)]

View(res.positive$pval.adj)


theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

names(res.positive)
dev.off()
# ENS ids seems to be required for the enrichment 
plot_enrichment_heatmap(res.positive)
ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)
dev.off()

plot_enrichment_heatmap(res.negative)
ggsave(paste0(outdir,'Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=120)


# Positive Factor 1: Hypoxia, oxygen dependent etc.
names(which(res.positive$pval.adj[,'Factor1']<0.0005))[1:10]
names(which(res.negative$pval.adj[,'Factor1']<0.05))

# Positive factor 2: citric acid cycle, meabolism, L13a-mediated, respratory
names(which(res.positive$pval.adj[,'Factor2']<0.00000005))[1:10]
names(which(res.negative$pval.adj[,'Factor2']<0.05))


# Positive factor 3: chromatin modifying enzymes, regulation of tp53
names(which(res.positive$pval.adj[,'Factor3']<0.05))[1:10]
names(which(res.negative$pval.adj[,'Factor3']<0.05))

#  
names(which(res.positive$pval.adj[,'Factor4']<0.05))[1:5]
names(which(res.negative$pval.adj[,'Factor4']<0.05))[1:5]






##### Make predictions 
## Here we select the top factors associated with the clinical variable of interest. 
## then we can choose the top variables of those factors 

## Prediction of clinical subgroups 
## Predict the EORTC.risk

suppressPackageStartupMessages(library(randomForest))

# Prepare data
df <- as.data.frame(get_factors(MOFAobject, factors=c(4))[[1]])
df


# Train the model for eortc.risk
df$EORTC.risk <- as.factor(MOFAobject@samples_metadata$EORTC.risk)
model.EORTC.risk <- randomForest(EORTC.risk ~ ., data=df, ntree=10)

# Do predictions
MOFAobject@samples_metadata$EORTC.risk.pred <- stats::predict(model.EORTC.risk, df)
MOFAobject@samples_metadata$EORTC.risk.pred

# Assess performance 
predicted<-as.factor(MOFAobject@samples_metadata$EORTC.risk.pred)
actual = MOFAobject@samples_metadata$EORTC.risk
confusion_mat = as.matrix(table(actual, predicted )) 
print(confusion_mat)

## Show "importance" of variables: higher value mean more important:
round(importance(model.EORTC.risk), 2)


# Prepare data
# Predict EORTC.risk with factor 1,2 only!
df <- as.data.frame(get_factors(MOFAobject, factors=c(3,4))[[1]])

# Train the model for IGHV
y_predict='EORTC.risk'
df$y <- as.factor(MOFAobject@samples_metadata[,y_predict])
model.y <- randomForest(y ~ .,data= df, ntree=10)
df$y <- NULL # important 
# Do predictions
MOFAobject@samples_metadata$y.pred <- stats::predict(model.y, df)
MOFAobject@samples_metadata$y.pred

# Assess performance 
predicted<-MOFAobject@samples_metadata$y.pred
actual <-as.factor(MOFAobject@samples_metadata[,y_predict])
confusion_mat = as.matrix(table(actual, predicted )) 
print(confusion_mat)
round(importance(model.y), 2)



### Plot predictions 






p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "EORTC.risk.pred",
                  shape_by = "EORTC.risk.pred_logical",
                  dot_size = 2.5,
                  show_missing = T
)

