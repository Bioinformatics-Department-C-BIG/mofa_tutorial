

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
rescale_option=FALSE
diff_var='NP2PTOT_LOG'
sel_group_cors = FALSE # Where to get corelations from
if (sel_group_cors){
  MOFAobjectPD_cors<-subset_groups(MOFAobjectPD, sel_group_cors)
}else{
  MOFAobjectPD_cors<-MOFAobjectPD
}

subset_groups

cors_both_clinical<-get_correlations(MOFAobjectPD_cors, all_diff_in_cors)
cors_pearson_pd_clinical = as.data.frame(cors_both_clinical[[2]]);  cors_all_pd_clinical = as.data.frame(cors_both_clinical[[1]])
## choose if the clu
all_fs_diff_all_time<-as.data.frame(cors_all_pd_clinical[,all_diff_in_cors]>(-log10(0.05)))
all_fs_diff = all_fs_diff_all_time
all_fs_diff[, c('NP2PTOT_LOG', 'updrs2_score_LOG_diff_V12')]

group_names(MOFAobject)

### Select group for plotting
sel_group=4
if (length(groups_names(MOFAobject))>1){
  MOFAobject_sel<-subset_groups(MOFAobject, sel_group)
  MOFAobjectPD_sel<-subset_groups(MOFAobjectPD, sel_group)

  met<-samples_metadata(MOFAobject_V08)
}else{
  MOFAobjectPD_sel = MOFAobjectPD
  MOFAobject_sel = MOFAobject
  met<-samples_metadata(MOFAobject)

}


all_fs_diff$NP2PTOT_LOG

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
        clusters_x=cluster_by_mofa_factors(MOFAobject=MOFAobjectPD_sel, centers=k_centers_m,
         factors=fact, rescale=rescale_option) # Cluster using kmeans
        clusters_x
       # if (MOFAobject==1){
        #  names(clusters_x)<-gsub('\\_.*', '', names(clusters_x))

       # }

        
        clust_ps<-clusters_x;clust_ps
        MOFAobjectPD_sel@samples_metadata[, xname] = as.factor(clust_ps)
        # TODO: save labels

        #write.csv()
        return(clust_ps)

}}
)
all_clusts_mofa$NP2PTOT_LOG
all_clusts_mofa_true<-all_clusts_mofa[!as.logical(lapply(all_clusts_mofa, is.null))]

all_clusts_mofa_true_t<-do.call(cbind,all_clusts_mofa_true)

dir.create(paste0(outdir,'/clustering/'))
all_clusts_file<-paste0(outdir,'/clustering/all_clusts_mofa.csv')

write.csv(all_clusts_mofa_true_t,all_clusts_file, row.names=TRUE,sep=','  )
  clusters_ids<-all_clusts_mofa[['NP2PTOT_LOG']]
  clusters_ids
# Attach clusters 
### Add all clusterings to mofa object 
diff_var = 'NP2PTOT_LOG'
for (diff_var in names(all_clusts_mofa)){
    print(diff_var)
    MOFAobjectPD_sel@samples_metadata[,paste0(diff_var, '_clust')]<-all_clusts_mofa[[diff_var]];#all_clusts_mofa[[diff_var]]
    
    
    #
    sm<-MOFAobject_sel@samples_metadata
    clusters_ids<-all_clusts_mofa[[diff_var]]
    clusters_ids
    rownames(clusters_ids )
    MOFAobject_sel@samples_metadata[,paste0(diff_var, '_clust')]<-clusters_ids[match(sm$PATNO_EVENT_ID,rownames(clusters_ids ) )]
    MOFAobject_sel@samples_metadata[(sm$INEXPAGE %in% c('INEXHC')),paste0(diff_var, '_clust')]<-'HC'
}


MOFAobject_sel@samples_metadata[, 'NP2PTOT_LOG_clust']



### Boxplots by cluster 
diff_variables=c('NP2PTOT_LOG')
diff_variables=c('NP2PTOT_LOG','NP2PTOT_diff_V13_V14', 'updrs2_score_LOG_diff_V12', 'NP2PTOT_V10')

diff_variables_to_p=c('NP2PTOT', 'scopa', 'updrs3_score')#, 'AGE' )
diff_variables_to_p=c('updrs2_score', 'scopa', 'updrs3_score')#, 'AGE' )
diff_variables_to_p=c('NP2PTOT_LOG', 'NP2PTOT','scopa','updrs3_score', 
                      'tremor','NP3BRADY', 'rigidity',   'rem', 'moca',
                      'upsit', 'VLTFRUIT', 'sft', 
                      'stai_state', 'stai_trait', 
                      'AGE_SCALED' )
                      


other_metrics<-t(cors_all_pd[all_fs_diff[,y],])
other_metrics[rowSums(other_metrics)>0,]
diff_variables_to_p=c('NP2PTOT', 'Neutrophil.Lymphocyte', 'AGE_SCALED', 'scopa', 
         'Neutrophils....', 'Lymphocytes....' ,   'sft', 'hvlt_immediaterecall',  # current 
         #baseline
         'NP2PTOT_BL',         
         # V10/12 -future scores 
         'NP2PTOT_V10',  'VLTFRUIT_V10', 'scopa_V10',
         'updrs3_score_LOG_V12',
          'hvlt_immediaterecall_V12'

         )

met$NP2PTOT_V10
 diff_variables_to_p=c( 
         'NP2PTOT',  'scopa', 'sft', 
          # for cognition#'moca',# 'hvlt_immediaterecall',  # current 
         #baseline
         'NP2PTOT_BL',         
         # V10/12 -future scores 
         'NP2PTOT_V10',  'VLTFRUIT_V10', 'scopa_V10', # 'moca_V12',
          #covars
         'Neutrophil.Lymphocyte', 'AGE_SCALED'
      #    'hvlt_immediaterecall_V12'

         )

                     #'Usable_Bases_SCALE' #, 'AGE' 
all_fs_diff[,'NP2PTOT_LOG']

MOFAobject

### Select group for plotting
met<-samples_metadata(MOFAobject_sel)
sel_group=4


y_clust="NP2PTOT_LOG"
y_clust="NP2PTOT_diff_V13_V14"
y_clust="updrs2_score_LOG_diff_V12"



diff_variables
facet_rows = 2
clust_name
sapply(diff_variables, function(y_clust){
  clust_name = paste0(y_clust, '_clust')
  ## check if there are clusters for this variable

  if (clust_name %in% colnames(met)){
    # for each cluster 


     # k centers might be different for each score 
   # k_centers <- max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)
    cluster_params<-paste0(clust_name ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors) )
    cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );
    print(cluster_params_dir)
    bn_all<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,  '.png')

    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,
    diff_variables_to_p, width=5+length(diff_variables_to_p)/facet_rows , bn=bn_all, facet_rows = facet_rows)


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

#k_centers <- max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)

cluster_params<-paste0(clust_name ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors) )
cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );


outfile_clusters<-paste0(cluster_params_dir, '/factor_plot_clusters_g' ,sel_group,  '.png')
color_by = 'NP2PTOT_clust'
outfile_clusters
# Plot clustering for scales 


p <- plot_factors(MOFAobjectPD_sel, 
             factors=which(all_fs_diff[,y]),
             color_by =color_by
#             shape_by = color_by
)
p
ggsave(outfile_clusters, width = 4, height = 4 )


## 2. Get the heatmap with the averages of all clusters 
## or show histograms of all clusters 
cors_all_pd[, 'NP2PTOT_LOG']

### Means by group 
library(dplyr)
sm$NP2PTOT_LOG_V10
diff_variables_to_p
diff_variables_to_p=c('NP2PTOT_LOG', 'NP2PTOT','scopa','updrs3_score', 
                      'tremor','NP3BRADY', 'rigidity',   'rem', 'moca',
                      'upsit', 'VLTFRUIT', 'sft', 
                      'stai_state', 'stai_trait', 
                      'AGE_SCALED', 'Neutrophil.Lymphocyte', 
                       'nfl_serum')
                      # nfl serum is the other way round why?
                      # Also check if the neutrophil: lymphocyte is higher in more severe patients 
                      # is cluster 3 with higher neutrophils more severe in the longterm or less?  


outfile_cl_heatmap<-paste0(cluster_params_dir, '/heatmap_means' ,  '.png')
diff_variables_to_p %in% colnames(samples_metadata(MOFAobjectPD_sel))
clust_name
col_data<-samples_metadata(MOFAobjectPD_sel)[c(diff_variables_to_p,clust_name, 'PATNO')]
colnames(col_data)
col_data$cluster<-col_data[, clust_name];
col_data[, clust_name]<-NULL
col_data$nfl_serum<-as.numeric(col_data$nfl_serum)

# Get the median per cluster and scale it for the ehatmap 

library(dplyr)
col_data_t<-tibble(col_data)
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
   labels_row= seq(1:k_centers_m),
  clustering_method = 'centroid',
  color = colorRampPalette(c("turquoise4", "white", "red"))(50))
dev.off()




#### TODO: tune with silhouette score ####


#### TODO: plot also clinical trajectories over time ####

source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' ))





## TODO: get heatmaps time v08 only 


### TODO: get network plots -
## TODO: add to the pipeline here the deseq groups or load them from the file!! 

  
if (FALSE){

selected_covars3 = selected_covars2

#corrplot::corrplot(stat[sig,], tl.col = "black", title="Pearson correlation coefficient")

#selected_covars_broad

selected_covars3<-selected_covars2[selected_covars2 %in% colnames(samples_metadata(MOFAobject))]
x1_all<-samples_metadata(MOFAobject)[selected_covars3][,1:10]

x2<-as.numeric(samples_metadata(MOFAobject)$cluster)
#cor_clust <- psych::corr.test(samples_metadata(MOFAobject)$td_pigd_old_on ,as.numeric(samples_metadata(MOFAobject)$cluster), method = "pearson", adjust = "BH")
samples_metadata(MOFAobject)[selected_covars2]
#install.packages('DescTools')
#apply(x1_all,2  ,MutInf, x2=x2, base=2)


#for (i in 1: length(selected_covars3)){
#  x1<-samples_metadata(MOFAobject)[selected_covars3][,i]
#  print(paste(selected_covars3[i], MutInf(x1, x2)))
#}



  ############## Create boxplots by group #### 
  col_data<-samples_metadata(MOFAobjectPD)[c(diff_variables_to_p,clust_name, 'PATNO')]
  col_data_melt<-reshape::melt(col_data, id=c('PATNO', clust_name))

  ggplot(col_data_melt)+ 
    geom_boxplot(aes_string(y='value', group=clust_name))+
    facet_wrap(~variable,scales =  'free')
    




  group_by(col_data, cluster) %>%
    summarise(across(everything(), list(~var(., na.rm=TRUE)))
              
    )

  
}




































































































































































































