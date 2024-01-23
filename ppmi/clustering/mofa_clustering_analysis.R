

#install.packages('DescTools')
library(DescTools)
library('factoextra')
# Using an existing trained model on simulated data
#file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#model <- load_model(file)

# Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 
# cluster samples here 



get_factors_for_metric<-function(diff_var){

  # get associated factors, remove the ones related to confounding
  fact <- which(all_fs_diff[,diff_var])

        if (remove_cell_factors){
          fact<-fact[!(fact %in% fact_neutro_pd)]

        }
       
        return(fact)
        }


#### Create the MOFA clusters with the same K ####
k_centers_m=3
remove_cell_factors = FALSE
#diff_var=y;  # diff_var='NP2PTOT_diff_V16'
rescale_option=TRUE
#diff_var='NP2PTOT_LOG'
#DIFF_VAR='moca'
DIFF_VAR='NP2PTOT_LOG'

sel_group_cors = FALSE # Where to get corelations from
if (sel_group_cors){
  MOFAobjectPD_cors<-subset_groups(MOFAobjectPD, sel_group_cors)
}else{
  MOFAobjectPD_cors<-MOFAobjectPD
}


cors_both_clinical<-get_correlations(MOFAobjectPD_cors, all_diff_in_cors)
cors_pearson_pd_clinical = as.data.frame(cors_both_clinical[[2]]);  cors_all_pd_clinical = as.data.frame(cors_both_clinical[[1]])
## choose if the clu
all_fs_diff_all_time<-as.data.frame(cors_all_pd_clinical[,all_diff_in_cors]>(-log10(0.05)))
all_fs_diff = all_fs_diff_all_time
all_fs_diff[, c('NP2PTOT_LOG', 'updrs2_score_LOG_diff_V12', 'moca')]

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
cors<-as.data.frame(cors)

cors$Neutrophil.Score;cors$Neutrophil.Lymphocyte; cors$Lymphocytes....
cors_all_pd$Neutrophil.Score;cors_all_pd$Neutrophil.Lymphocyte
# Remove factors that have to do with site, plate, 
c((which(cors$`Lymphocytes....`>-log10(0.05))), which(cors$`Neutrophils....`>-log10(0.05)), which(cors$`SITE`>-log10(0.05)), which(cors$`Plate`>-log10(0.05)) )


fact_neutro_pd<-c(  (which(cors_all_pd$`Lymphocytes....`>-log10(0.05))), which(cors_all_pd$`Neutrophils....`>-log10(0.05)),  which(cors_all_pd$Plate>-log10(0.05)), 
                      which(cors_all_pd$RIN>-log10(0.05)), which(cors_all_pd$`Uniquely.mapped....`>-log10(0.05)),
                       which(cors_all_pd$Usable_bases_SCALE>-log10(0.05))
                      #which(cors$`Uniquely.mapped....`>-log10(0.05))
                        )

fact_neutro_pd

fact_neutro_pd<-c(  (which(cors$`Lymphocytes....`>-log10(0.05))), which(cors$`Neutrophils....`>-log10(0.05)),  which(cors$Plate>-log10(0.05)), 
                      which(cors$RIN>-log10(0.05)), which(cors$Usable_bases_SCALE>-log10(0.05)))
                      #which(cors$`Uniquely.mapped....`>-log10(0.05))
                        
fact_neutro_pd <-c(fact_neutro_pd, 8,2)
fact_neutro_pd<-c(2)
DIFF_VAR
    
get_factors_for_metric(DIFF_VAR)


all_clusts_mofa <- sapply(colnames(all_fs_diff),function(diff_var){
  #'
  #' @param
  #' 
  #' 

        fact <- get_factors_for_metric(diff_var)
    


        print(paste(diff_var,paste(fact, collapse=', ')))
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





for (diff_var in names(all_clusts_mofa)){
    MOFAobjectPD_sel@samples_metadata[,paste0(diff_var, '_clust')]<-all_clusts_mofa[[diff_var]];#all_clusts_mofa[[diff_var]]
    sm<-MOFAobject_sel@samples_metadata
    clusters_ids<-all_clusts_mofa[[diff_var]]
    clusters_ids
    rownames(clusters_ids )
    MOFAobject_sel@samples_metadata[,paste0(diff_var, '_clust')]<-clusters_ids[match(sm$PATNO_EVENT_ID,rownames(clusters_ids ) )]
    MOFAobject_sel@samples_metadata[(sm$INEXPAGE %in% c('INEXHC')),paste0(diff_var, '_clust')]<-'HC'
}


MOFAobject_sel@samples_metadata[, 'NP2PTOT_LOG_clust']



### Boxplots by cluster 
## Can produce for multiple metrics 
diff_variables=c('NP2PTOT_LOG')
diff_variables=c('scopa', 'moca', 'NP2PTOT_LOG','NP2PTOT_diff_V13_V14', 'updrs2_score_LOG_diff_V12', 'NP2PTOT_LOG_V10', 'moca')

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


# 1. get the factors 
cors_all_pd<-as.data.frame(cors_all_pd)
get_covariates_cells<-function(y_clust, thresh=0){
  # get the cell types that corelated with the metric
  #' 
  #' @thresh: -log10pvalue
  fact=get_factors_for_metric(y_clust)
  colnames(estimations)[!colnames(estimations) %in% colnames(cors_all_pd)]
  
  cell_types<-c(colnames(estimations), measured_cells)[c(colnames(estimations), measured_cells) %in% colnames(cors_all_pd) ]
  
  c(colnames(estimations), measured_cells)

  clust_variates<-cors_all_pd[fact, cell_types]
  clust_variates
  variates_to_p<-names(which(colSums(clust_variates)>thresh))
  variates_to_p
  return(variates_to_p)
}

get_covariates_cells('NP2PTOT_LOG')


paste0(DIFF_VAR,'_BL')

 diff_variables_to_p=c( 
         'NP2PTOT',  'scopa', 'sft', 
          # for cognition#'moca',# 'hvlt_immediaterecall',  # current 
         #baseline
        # 'NP2PTOT_BL',     
        paste0(DIFF_VAR,'_BL'),     
         # V10/12 -future scores 
        # 'VLTFRUIT_V10', 'scopa_V10', # 'moca_V12',
#         paste0(DIFF_VAR,'_V10'), paste0(DIFF_VAR,'_V12'),  paste0('sft_V12'), 
        'NP2PTOT_V10', 'NP2PTOT_V12','sft_V12',
        # 'MCATOT_V14',
          #covars
         'Neutrophil.Lymphocyte', 'AGE_SCALED', 'moca',
         'duration'
      #    'hvlt_immediaterecall_V12'

         )

### Select group for plotting
met<-samples_metadata(MOFAobject_sel)
sel_group=4


y_clust="NP2PTOT_LOG"
y_clust=DIFF_VAR
diff_variables
clust_name
facet_rows = 2


## BOXPLOTS 

sapply(diff_variables, function(y_clust){
  clust_name = paste0(y_clust, '_clust')
  ## check if there are clusters for this variable

  if (clust_name %in% colnames(met)){
    # for each cluster create boxplot 
    variates_to_p<-get_covariates_cells(y_clust, thresh=8)
    fact=get_factors_for_metric(y_clust)
    fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))
#  'rcf_',as.numeric(remove_cell_factors )
    cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors))
    cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );

    bn_all<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,  '.png')

    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,
    c(diff_variables_to_p, variates_to_p), width=8+length(c(diff_variables_to_p,variates_to_p))/facet_rows , bn=bn_all, facet_rows = facet_rows, 
    text=paste('removed',paste0( unique(fact_neutro_pd),  collapse=', ')))


    }else{
      print(paste0(clust_name, '_clust does not exist'))
    }
  })


### Variables in this script 
# 1. Plot the clusters on the factor plot 
#' @param all_fs_diff # table of clinical scores and factors: which factors are sign with which score
DIFF_VAR
y <- DIFF_VAR# cluster metric 
DIFF_VAR
color_by=paste0(y, '_clust')
clust_metric<-y
y
# Settings for each clustering 

met<-met[!is.na(met[, clust_metric]),]
clust_name<-paste0(clust_metric, '_clust')
print(paste('Using subset of  ', dim(met)[1], ' patients'))
freqs<-paste0('n=', paste0(table(met[, clust_name]), collapse = ', '))

#k_centers <- max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)

fact=get_factors_for_metric(y_clust)

fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))

cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors) )
cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );


outfile_clusters<-paste0(cluster_params_dir, '/factor_plot_clusters_g' ,sel_group, y, '_', color_by, '.png')
outfile_clusters
# Plot clustering for scales 


p <- MOFA2::plot_factors(MOFAobjectPD_sel, 
             factors=which(all_fs_diff[,y]),
             color_by =color_by
            # alpha=0.7
#             shape_by = color_by
)
p
ggsave(outfile_clusters, width = 4, height = 4 )
dev.off()

## 2. Get the heatmap with the averages of all clusters 
## or show histograms of all clusters 


### Means by group 
library(dplyr)



diff_variables_to_p=c(diff_variables_to_p, 'nfl_serum')
diff_variables_to_p
#diff_variables_to_p=c('NP2PTOT_LOG', 'NP2PTOT','scopa','updrs3_score', 
#                      'tremor','NP3BRADY', 'rigidity',   'rem', 'moca',
#                      'upsit', 'VLTFRUIT', 'sft', 
#                      'stai_state', 'stai_trait', 
#                      'AGE_SCALED', 'Neutrophil.Lymphocyte', 
#                       'nfl_serum')
                      # nfl serum is the other way round why?
                      # Also check if the neutrophil: lymphocyte is higher in more severe patients 
                      # is cluster 3 with higher neutrophils more severe in the longterm or less?  


outfile_cl_heatmap<-paste0(cluster_params_dir, '/heatmap_means' ,  '.png')
diff_variables_to_p %in% colnames(samples_metadata(MOFAobjectPD_sel))
clust_name=paste0(DIFF_VAR, '_clust')

col_data<-samples_metadata(MOFAobjectPD_sel)[c(diff_variables_to_p,clust_name, 'PATNO')]
colnames(col_data)
col_data$cluster<-col_data[, clust_name];
col_data[, clust_name]<-NULL
col_data$nfl_serum<-as.numeric(col_data$nfl_serum)

# Get the median per cluster and scale it for the ehatmap 

library(dplyr)
col_data_t<-tibble(col_data)
#means_by_cluster %>% group_indices()
colnames(col_data)
means_by_cluster <- col_data_t %>% 
      dplyr::select(-PATNO) %>%
      group_by(cluster) %>%
      mutate_if(is.character, as.numeric) %>%
  summarise_each(funs(median(., na.rm = TRUE)))%>%
  scale() %>% as.data.frame() %>% dplyr::select(-cluster) %>% 
  as.data.frame()
rownames(means_by_cluster)<-means_by_cluster$cluster


# TODO: adjust clustering method to get the correct roder
# reverse the moca scores because the severe is low
means_by_cluster_na<-means_by_cluster[,colSums(!is.na(means_by_cluster))>0]
dim(means_by_cluster_na)
mc_scores<-grep('MC|moca|sft',colnames(means_by_cluster_na))
means_by_cluster_na[,mc_scores]=-means_by_cluster_na[,mc_scores]
colnames(means_by_cluster_na)[mc_scores] = paste0(colnames(means_by_cluster_na)[mc_scores] , '_reverse')
outfile_cl_heatmap
graphics.off()
jpeg(outfile_cl_heatmap, width=7, height=3, units='in', res=200)

  pheatmap(means_by_cluster_na, 
   labels_row= seq(1:k_centers_m),
  clustering_method = 'centroid',
  color = colorRampPalette(c("turquoise4", "white", "red"))(50))
dev.off()
graphics.off()




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

sm$duration




