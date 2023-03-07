
#install.packages('psych')

#install.packages("remotes")
#remotes::install_github("bioFAM/MOFAdata")
##### 
### this one depends on mofa application to inherit: 
### 1. MOFAobject, 2. factors, 
# 3. clinical variables
# 4. outdirs 
# 5. csf/plasma/untargeted flags


#source('enrichment.R')
library(ggplot2)




print(outdir)
dev.off()

jpeg(paste0(outdir, 'factor_cor','.jpeg'))
plot_factor_cor(MOFAobject)
dev.off()


calculate_variance_explained(MOFAobject)
plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained','.png'), width = 4, height=4, dpi=100)



MOFAobject@samples_metadata$NP1ANXS<-as.factor(MOFAobject@samples_metadata$NP1ANXS)
MOFAobject@samples_metadata$Grade
colnames(MOFAobject@samples_metadata)[20:35]


stats<-apply(MOFAobject@samples_metadata, 2,table )
sapply(stats,var)>0

non_na_vars<-which(!is.na(sapply(stats,mean)) & sapply(stats,var)>0 )

NROW(non_na_vars)
#### Covariance of factors with metadata 
correlate_factors_with_covariates(MOFAobject,
                                  covariates = c('INFODT.x','INFODT.y' ,'NUPSOURC',"NP1ANXS", "NP1DPRS", 'NP1HALL', 'NP1APAT', 'NP1DDS', 'NP1RTOT'), 
                                  plot = "log_pval",
                                  
)

correlate_factors_with_covariates(MOFAobject,
                                  covariates = c("NP1ANXS", "NP1DPRS", 'NP1HALL', 'NP1APAT', 'NP1DDS', 'NP1RTOT'), 
                                  plot = "log_pval"
                                 
                                  
)
dev.off()



## plot only factors that correlate? 
cors<-correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(non_na_vars), 
                                  plot = "log_pval", 
                                  return_data = TRUE
                                  
)
ids_to_plot<-which(apply(cors, 2, sum)>0)



jpeg(paste0(outdir, 'factors_covariates_all','.jpeg'), width = 2000, height=700, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(non_na_vars), 
                                  plot = "log_pval"
                                  
)
dev.off()


jpeg(paste0(outdir, 'factors_covariates_only_nonzero','.jpeg'), width = 2000, height=700, res=150)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = names(non_na_vars)[ids_to_plot], 
                                  plot = "log_pval"
                                  
)
dev.off()

### filter only the ones that are correlated 

covariate_corelations<-correlate_factors_with_covariates(MOFAobject,
                                                         covariates = colnames(MOFAobject@samples_metadata)[c(6:12,45:70, 90:124)], 
                                                         plot = "log_pval",
                                                         return_data = TRUE
)

covariate_corelations<-correlate_factors_with_covariates(MOFAobject,
                                  covariates = colnames(MOFAobject@samples_metadata)[c(6:12,45:70, 90:124)], 
                                  plot = "log_pval",
                                  return_data = TRUE
)
write.csv(covariate_corelations, paste0(outdir, '/covariate_corelations.csv'))

jpeg(paste0(outdir, 'factors_covariates','.jpeg'), width = 2000, height=800, res = 200)
correlate_factors_with_covariates(MOFAobject,
                                  covariates = colnames(MOFAobject@samples_metadata)[c(6:12,45:65, 90:124)], 
                                  plot = "log_pval"
)
dev.off()

view='proteomics'; factor=6

vps=length(MOFAobject@dimensions$D)
fps= as.numeric(MOFAobject@dimensions$K)
fps
views<-names(MOFAobject@dimensions$D)
views


##### WRITE ALL weights for each factor in one file 

for (i in seq(1,vps)){
  view=views[i]
  all_weights<-MOFA2::get_weights(MOFAobject,views = view, 
                                  as.data.frame =TRUE)  
  # threshold each? 
  T=0.3
  all_weights_filt<-all_weights[abs(all_weights$value)>T,]
  dim(all_weights_filt)
  write.table(all_weights_filt,paste0(outdir, 'top_weights_vals_by_view_', view, '_T_', T, '.txt'), sep = '\t')
  ens_ids<-gsub('\\..*', '', all_weights_filt$feature)
  write.csv(ens_ids,paste0(outdir, 'top_weights_vals_by_view_CLUEGO_', view, '_T_', T, '.txt'),
            row.names = FALSE, quote=FALSE)
  print(view)
  print(dim(all_weights_filt))
  print(dim(all_weights))
  
  }




for (i in seq(1,vps)){
  for (ii in seq(1,fps)){
    view=views[i]
    factor=ii
    print(view, factor)
    all_weights<-MOFA2::get_weights(MOFAobject,views = view, factors=factor, 
                            as.data.frame =TRUE)
    
    ### get the top highly weighted variables - absolute value
    top<-all_weights[order(abs(all_weights$value), decreasing = TRUE),]
    write.table(top,paste0(outdir, 'top_weights_vals',factor,'_', view,'.txt'), sep = '\t')
    

  }
  }

  
  plot_variance_explained(MOFAobject, plot_total = T)[[2]]
  ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)
#install.packages('psych')

  
  
#### Get all weights and put it in one file
  #1. Collate in a list  - order significant
  #2. Stack the lists - 
  

  
  
  
  
  
  
  
  
  
  
  
plot_factors(MOFAobject, 
             factors = c(3,8), 
             dot_size = 2.5, 
             color_by = 'NP3BRADY'
             
)


##### Plot molecular signatures in the input data



plot_weights(MOFAobject,
             view = "miRNA",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)





### Age, gender , stage does not discriminate factors
#Conclusion here:# factor 1 correlates with grade

p<-plot_factor(MOFAobject, 
            factors = 3, 
            color_by = "NP3BRADY",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)
p
plot_factor(MOFAobject, 
            factors = 4, 
            color_by = "NP3BRADY",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "NP3BRADY",
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)



# Factor 2 associates with proteomic Subtype 

color_by<-'NHY';fs<-c(4,8)
color_by<-'NP3TOT';fs<-c(1,4)


plot_factor(MOFAobject, 
            factors = fs, 
            color_by = color_by,
            add_violin = TRUE,
            dodge = TRUE,
            show_missing = FALSE
)
fss<-paste(fs,sep='_',collapse='-')
ggsave(paste0(outdir, 'plot_factor_variate_violin',fss,color_by,'.png'), width = 4, height=4, dpi=100)

#### plot 2 factors 
##### TODO: Plot only significant covariates  here

color_by<-'NHY'
fs<-c(4,8)

color_by<-'NHY';fs<-c(1,8)
color_by<-'NP3TOT';fs<-c(1,8)

plot_factors(MOFAobject, 
            factors = fs, 
            color_by = color_by,
            show_missing = FALSE
)
fss<-paste(fs,sep='_',collapse='-')
FNAME<-paste0(outdir, 'plot_factors_variate_2D',fss,color_by,'.png')
FNAME
color_by
ggsave(FNAME, width = 4, height=4, dpi=100)



##### plot weights 

v_set=c()
v_set=c()

view='miRNA'
factor=8
plot_top_weights(MOFAobject,
                 view = view,
                 factor = factor,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
ggsave(paste0(outdir, 'top_weights_',factor, view,'_','.png'), width =3 , height=4, dpi=100)


fps=8
for (i in 1:vps){
  for (ii in 1:fps){
    print(c(i,ii))
    plot_top_weights(MOFAobject,
                     view = views[i],
                     factor = ii,
                     nfeatures = 10,     # Top number of features to highlight
                     scale = T           # Scale weights from -1 to 1
    )
    ggsave(paste0(outdir, 'top_weights_', ii,'_',vps[i],'.png'), width = , height=4, dpi=100)
    
    plot_weights(MOFAobject, 
                 view = views[i], 
                 factor = ii, 
                 nfeatures = 10
    )
    ggsave(paste0(outdir, 'all_weights_', ii,'_',vps[i],'.png'), width = 4, height=4, dpi=100)
    
    
    
    ###### Heatmaps 
    nfs=40
    print('heatmap')
    jpeg(paste0(outdir, 'heatmap_',ii,'_',views[i], 'nfs_', nfs, '.jpeg'), res=150,height=20*nfs, width=20*nfs)
    fps[ii]=1
    # Plot heatmaps for each factor only for miRNA 
    p<-plot_data_heatmap(MOFAobject, 
                         view = views[i], 
                         factor =  ii,  
                         features = nfs,
                         denoise = TRUE,
                         cluster_rows = TRUE, cluster_cols = FALSE,
                         show_rownames = TRUE, show_colnames = FALSE,
                         scale = "row",
                         fontsize_number = 5
                         
                         
    )
    dev.off()
    
    jpeg(paste0(outdir, 'heatmap_',ii,'_',views[i], 'nfs_', nfs, '.jpeg'), height=20*nfs, width=20*nfs)
    
    p<-plot_data_heatmap(MOFAobject, 
                         view = views[i], 
                         factor =  ii,  
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




# plot heatmaps by view and 



plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 1,  
                  features = 200,
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
                  factor = 2,  
                  features = 30,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

plot_data_heatmap(MOFAobject, 
                  view = "proteomics",
                  factor = 3,  
                  features = 100,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


### So multi omics factors are more related to Stage than to subtype!!
##How?? plot on the factors and color by stage

p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "NP3BRADY",
                  shape_by = "NP3BRADY",
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

#BiocManager::install('AnnotationHub')

source('enrichment.R')
  
#library(AnnotationHub)
#ah = AnnotationHub()


library('MOFAdata')
utils::data(reactomeGS)

head((reactomeGS))



subcategory<- 'CP:KEGG'

subcategory<- 'GO:MF'
subcategory<- 'GO:BP'

gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
gs<-as.matrix(read.csv(gs_file, header=1, row.names=1))


features_names(MOFAobject)$RNA<-sapply(features_names(MOFAobject)$RNA, 
       function(x) {stringr::str_remove(x, '\\..*')}
)


# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "RNA",
                               sign = "positive"
)



# GSEA on negative weights, with default options
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = gs, 
                               view = "RNA",
                               sign = "negative"
)



res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = gs, 
                               view = "RNA",
                               sign = "positive"
)



  
  
  
# change to negative and positive

extract_order_significant<-function(x) {
  # extracy most significant and order 
  T=0.005
  sign<-x[x<T]
  sign2<-sign[order(sign)]
  print(sign2)
}

enrichment_list=all_fs_enrichment
stack_list<-function(i,enrichment_list) {
  
  # Take a list of dataframes and stack them 
  # Add a column called factor which extracts the list counter 
  
  x=enrichment_list[[i]]
  if (length(x)){
    tmp<-as.data.frame(x)
    tmp$path<-rownames(tmp)
    colnames(tmp)<-'pvals'
    
    rownames(tmp)<-NULL
    f<-names(all_fs_enrichment)[[i]] # EXTRACT the counter to assign factor value in new column
    tmp$factor=f
    return(tmp)}}



results_enrich<-res.positive$pval.adj
all_fs_enrichment<-apply(results_enrich, 2 , extract_order_significant)

all_fs_unlisted<-sapply(seq(1:length(all_fs_enrichment)), stack_list, enrichment_list=all_fs_enrichment)
all_fs_merged<-do.call(rbind, all_fs_unlisted )
dim(all_fs_merged)




write.csv(unlist(all_fs_enrichment, use.names = TRUE),paste0(outdir, gsub('\\:', '_', subcategory), '_enrichment_positive_pvals', out_params, '.csv' ))
write.csv(all_fs_merged,paste0(outdir, gsub('\\:', '_', subcategory), '_enrichment_positive_pvals_no_f_' ,  out_params, '.csv' ))



results_enrich<-res.negative$pval.adj
all_fs_enrichment<-apply(results_enrich, 2 , extract_order_significant)
all_fs_unlisted<-sapply(seq(1:length(all_fs_enrichment)), stack_list, enrichment_list=all_fs_enrichment)
all_fs_merged<-do.call(rbind, all_fs_unlisted )
dim(all_fs_merged)

write.csv(unlist(all_fs_enrichment, use.names = TRUE),paste0(outdir,gsub('\\:', '_', subcategory), '_enrichment_negative_pvals', out_params, '.csv' ))
write.csv(all_fs_merged,paste0(outdir,gsub('\\:', '_', subcategory), '_enrichment_negative_pvals_no_f_', out_params, '.csv' ))





# Make enrichment plots for all factors 
# threshold on p value to zoom in 
jpeg(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.positive, 
                        alpha=0.5, cap=0.0005)
dev.off()

plot_enrichment_heatmap(res.positive$sigPathways, 
                        alpha=0.5, cap=0.0005)

#ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)


jpeg(paste0(outdir,'Enrichment_heatmap_negative','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.negative, 
                        alpha=0.5, cap=0.0005)
dev.off()
#ggsave(paste0(outdir,'Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=120)


F3<-res.positive$pval.adj[,'Factor3']
SIG<-F3[F3<0.05]
SIG[order(SIG)][1:20]

F3<-res.negative$pval.adj[,'Factor6']
SIG<-F3[F3<0.05]
SIG[order(SIG)][1:10]



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
y_predict='NHY'
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

