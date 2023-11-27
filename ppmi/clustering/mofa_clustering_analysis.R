

#install.packages('DescTools')
library(DescTools)
library('factoextra')
# Using an existing trained model on simulated data
#file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#model <- load_model(file)

# Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 
# cluster samples here 

#### Create the MOFA clusters with the same K ####
k_centers_m=3
#diff_var=y;  # diff_var='NP2PTOT_diff_V16'
all_fs_diff$NP2PTOT_LOG
rescale_option=TRUE
all_clusts_mofa <- sapply(colnames(all_fs_diff),function(diff_var){
  #'
  #' @param
  #' 
  #' 

        fact <- which(all_fs_diff[,diff_var])
        xname = paste0(diff_var, '_clust')
      if (length(fact) > 0){
        set.seed(42)
        set.seed(60)
        clusters_x=cluster_by_mofa_factors(MOFAobject=MOFAobjectPD, centers=k_centers_m,
         factors=fact, rescale=rescale_option) # Cluster using kmeans
        clust_ps<-clusters_x$cluster
        MOFAobjectPD@samples_metadata[, xname] = as.factor(clust_ps)
        MOFAobjectPD@samples_metadata[, xname] = as.factor(clust_ps)

        return(clust_ps)
}}
)

# Attach clusters 
### Add all clusterings to mofa object 
for (diff_var in names(all_clusts_mofa)){
    MOFAobjectPD@samples_metadata[,paste0(diff_var, '_clust')]<-all_clusts_mofa[[diff_var]]
    sm<-MOFAobject@samples_metadata
    clusters_ids<-all_clusts_mofa[[diff_var]]
    MOFAobject@samples_metadata[,paste0(diff_var, '_clust')]<-clusters_ids[match(sm$PATNO,names(clusters_ids ) )]
    MOFAobject@samples_metadata[(sm$INEXPAGE %in% c('INEXHC')),paste0(diff_var, '_clust')]<-'HC'
}
MOFAobject@samples_metadata$NP2PTOT_LOG_clust






### Boxplots by cluster 
diff_variables=c('NP2PTOT_LOG')
diff_variables_to_p=c('NP2PTOT', 'scopa', 'updrs3_score')#, 'AGE' )
diff_variables_to_p=c('updrs2_score', 'scopa', 'updrs3_score')#, 'AGE' )
diff_variables_to_p=c('NP2PTOT_LOG','scopa','updrs3_score', 
                      'tremor','NP3BRADY', 'rigidity',   'rem', 'moca',
                      'upsit', 'VLTFRUIT', 'sft', 
                      'stai_state', 'stai_trait', 
                      'AGE_SCALED' )
                     #'Usable_Bases_SCALE' #, 'AGE' 
all_fs_diff[,'NP2PTOT_LOG']

met<-samples_metadata(MOFAobject)
met$updrs2_score_LOG
y_clust="NP2PTOT_LOG"


sapply(diff_variables, function(y_clust){
  clust_name = paste0(y_clust, '_clust')
  ## check if there are clusters for this variable
  print(table(samples_metadata(MOFAobject)[, clust_name]))
  
  if (clust_name %in% colnames(met)){


     # k centers might be different for each score 
    k_centers<-max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)
    cluster_params<-paste0(clust_name ,'/', k_centers,'/',rescale_option)
    cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );
    cluster_params_dir

    bn<-paste0(cluster_params_dir , y,  '.png') # need to fix to add y in every iteration in function
     #   sapply(diff_variables_to_p_all, boxplot_by_cluster, met=met, clust_name=clust_name, bn=bn)

    bn_all<-paste0(cluster_params_dir, '/all_vars' ,  '.png')
    bn_all

    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,diff_variables_to_p, width=15, bn=bn_all)
    }
  })



### Variables in this script 
# 1. Plot the clusters on the factor plot 
#' @param all_fs_diff # table of clinical scores and factors: which factors are sign with which score
y <- 'NP2PTOT_LOG' # cluster metric 
color_by=paste0(y, '_clust')
clust_metric<-y
# Settings for each clustering 

met<-met[!is.na(met[, clust_metric]),]
clust_name<-paste0(clust_metric, '_clust')
print(paste('Using subset of  ', dim(met)[1], ' patients'))
freqs<-paste0('n=', paste0(table(met[, clust_name]), collapse = ', '))

k_centers <- max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)
cluster_params_dir <- paste0(outdir,'/clustering/',clust_name ,'/', k_centers,'/',rescale_option );
cluster_params_dir
outfile_clusters<-paste0(cluster_params_dir, '/factor_plot_clusters' ,  '.png')
color_by = 'NP2PTOT_clust'

MOFAobjectPD@samples_metadata$NP2PTOT
# Plot clustering for scales 
p <- plot_factors(MOFAobjectPD, 
             factors=which(all_fs_diff[,y]),
             color_by =color_by
#             shape_by = color_by
)
p
ggsave(outfile_clusters, width = 4, height = 4 )


## 2. Get the heatmap with the averages of all clusters 
## or show histograms of all clusters 


### Means by group 
library(dplyr)
diff_variables_to_p
diff_variables_to_p=c('NP2PTOT_LOG','scopa','updrs3_score', 
                      'tremor','NP3BRADY', 'rigidity',   'rem', 'moca',
                      'upsit', 'VLTFRUIT', 'sft', 
                      'stai_state', 'stai_trait', 
                      'AGE_SCALED' )
                      # nfl serum is the other way round why? 


outfile_cl_heatmap<-paste0(cluster_params_dir, '/heatmap_means' ,  '.png')
diff_variables_to_p %in% colnames(samples_metadata(MOFAobjectPD))
clust_name
col_data<-samples_metadata(MOFAobjectPD)[c(diff_variables_to_p,clust_name, 'PATNO')]

MOFAobjectPD@samples_metadata[, clust_name]
col_data$rem
col_data$cluster<-col_data[, clust_name];col_data[, clust_name]<-NULL

# Get the median per cluster and scale it for the ehatmap 

#means_by_cluster %>% group_indices()
means_by_cluster <- col_data %>% 
      dplyr::select(-PATNO) %>%
      group_by(cluster) %>%
      mutate_if(is.character, as.numeric) %>%
  summarise_each(funs(median(., na.rm = TRUE)))%>%
  scale() %>% as.data.frame() %>% select(-cluster) %>% 
  as.data.frame()
rownames(means_by_cluster)<-means_by_cluster$cluster


# TODO: adjust clustering method to get the correct roder
means_by_cluster_na<-means_by_cluster[,colSums(!is.na(means_by_cluster))>0]
dim(means_by_cluster_na)
means_by_cluster_na$moca_rev<- (-means_by_cluster_na$moca) # reverse moca
means_by_cluster_na$moca<-NULL
graphics.off()
jpeg(outfile_cl_heatmap, width=7, height=3, units='in', res=200)
  pheatmap(means_by_cluster_na, 
   labels_row= c(1,2,3),
  clustering_method = 'centroid',
  color = colorRampPalette(c("turquoise4", "white", "red"))(50))
dev.off()




#### TODO: tune with silhouette score ####


#### TODO: plot also clinical trajectories over time ####

source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' ))





## TODO: get heatmaps time v08 only 


### TODO: get network plots -
## TODO: add to the pipeline here the deseq groups or load them from the file!! 

  



corrplot::corrplot(stat[sig,], tl.col = "black", title="Pearson correlation coefficient")




x1_all<-samples_metadata(MOFAobject)[selected_covars3][,1:10]

x2<-as.numeric(samples_metadata(MOFAobject)$cluster)
cor_clust <- psych::corr.test(samples_metadata(MOFAobject)[selected_covars2]$td_pigd_old_on ,as.numeric(samples_metadata(MOFAobject)$cluster), method = "pearson", adjust = "BH")
samples_metadata(MOFAobject)[selected_covars2]
#install.packages('DescTools')
apply(x1_all,2  ,MutInf, x2=x2, base=2)


for (i in 1: length(selected_covars3)){
  x1<-samples_metadata(MOFAobject)[selected_covars3][,i]
  print(paste(selected_covars3[i], MutInf(x1, x2)))
}




############## Create boxplots by group #### 
col_data<-samples_metadata(MOFAobjectPD)[c(diff_variables_to_p,clust_name, 'PATNO')]
col_data_melt<-reshape::melt(col_data, id=c('PATNO', clust_name))

ggplot(col_data_melt)+ 
  geom_boxplot(aes_string(y='value', group=clust_name))+
  facet_wrap(~variable,scales =  'free')
  




group_by(col_data, cluster) %>%
  summarise(across(everything(), list(~var(., na.rm=TRUE)))
            
  )

 













































































