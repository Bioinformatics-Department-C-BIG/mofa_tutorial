
#### Outline ####
# 1. Boxplots
# 2. 

# 1. get the factors 


cors_all_pd<-as.data.frame(cors_all_pd)

# Find out the top correlates
df1<-cors_all_pd['Factor2',]%>%t() %>% as.data.frame()
df1 %>% dplyr::arrange(Factor2 )



paste0(DIFF_VAR,'_BL')

 diff_variables_to_p=c( 
         'NP2PTOT',
         'NP3TOT',  'scopa', 'sft',  'moca',
         'Neutrophil.Lymphocyte', 'AGE',
         'duration', 
         # f23
         # added V12 because it is significant ? 
         'RBD_TOT','rem_V12'
                  )



### Select group for plotting
met<-samples_metadata(MOFAobject_sel)
sel_group=4


y_clust="NP2PTOT_LOG"
y_clust=DIFF_VAR
clust_vars<-c('NP2PTOT_LOG', 'moca', 'NP3TOT_LOG', 'scopa', 'updrs3_score_on')

get_factors_for_metric('updrs3_score_on')
facet_rows = 2

#write_vars_output(MOFAobject, vars_by_factor, factors=fact)

#### Boxplots ####

sapply(clust_vars, function(y_clust){
  clust_name = paste0(y_clust, '_clust')
  ## check if there are clusters for this variable

  #clust_name %in% colnames(met)
  if (clust_name %in% colnames(met)){
    # for each cluster create boxplot 
    variates_to_p<-get_covariates_cells(y_clust, thresh=8)
    variates_to_p = c('B.Naive', 'B.Memory', 'Basophils.LD', 'Monocytes.C', 'Neutrophils.LD', 'T.CD4.Memory', 'T.CD4.Naive', 'T.CD8.Memory', 
    'T.CD8.Naive', 'Neutrophils....'  )
    variates_to_p<-variates_to_p[variates_to_p%in% colnames(met)]

    fact=get_factors_for_metric(y_clust)
    fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))

    cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors))
    cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );

    bn_all_fname<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,  '.png')


    diff_variables_to_p<-diff_variables_to_p[diff_variables_to_p %in% colnames(met)] 


    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,  c(diff_variables_to_p), width=8+length(c(diff_variables_to_p))/facet_rows, 
    height=1+1.5*facet_rows, bn=bn_all_fname, facet_rows = 1, 
    text='')


    bn_all_fname<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,'cell_types',  '.png')
    variates_to_p = colnames(estimations)[!colnames(estimations) %in% c('T.gd.Vd2',"Plasmablasts" , 'mDCs', 'Monocytes.NC.I' )] # TODO: check for zero variance to exclude 
    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,  c(variates_to_p), width=8+length(c(variates_to_p))/facet_rows , bn=bn_all_fname, facet_rows = facet_rows, 
    text='')


    }else{
      print(paste0(clust_name, '_clust does not exist'))
    }
  })



# 2. Plot the clusters on the factor plot  ####
#' @param all_fs_diff # table of clinical scores and factors: which factors are sign with which score
DIFF_VAR='NP2PTOT_LOG'
y <- DIFF_VAR# cluster metric 

color_by=paste0(y, '_clust')
clust_metric<-y

# Settings for each clustering 

met<-met[!is.na(met[, clust_metric]),]
clust_name<-paste0(clust_metric, '_clust')
freqs<-paste0('n=', paste0(table(met[, clust_name]), collapse = ', '))

#k_centers <- max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)

fact=get_factors_for_metric(y_clust)

fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))

cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors) )
cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );


outfile_clusters<-paste0(cluster_params_dir, '/factor_plot_clusters_g' ,sel_group, y, '_', color_by, '.png')

# Plot clustering for scales 


p <- MOFA2::plot_factors(MOFAobjectPD_sel, 
             factors=which(all_fs_diff[,y]),
             color_by =color_by
            # alpha=0.7
#             shape_by = color_by
)
p
ggsave(outfile_clusters, width = 4, height = 4 )




## 3. Get the heatmap with the averages of all clusters  ####
## or show histograms of all clusters 
### Means by group 
library(dplyr)



diff_variables_to_p=c(diff_variables_to_p, 'nfl_serum')


outfile_cl_heatmap<-paste0(cluster_params_dir, '/heatmap_means' ,  '.png')
diff_variables_to_p %in% colnames(samples_metadata(MOFAobjectPD_sel))


clust_name=paste0(DIFF_VAR, '_clust')
diff_variables_to_p
col_data<-samples_metadata(MOFAobjectPD_sel)[c(diff_variables_to_p,clust_name, 'PATNO')]
col_data$nfl_serum<-as.numeric(col_data$nfl_serum)
col_data$cluster<-col_data[, clust_name]; col_data[, clust_name]<-NULL




plot_heatmap_median_by_cluster<-function(col_data){
      #
    #means_by_cluster %>% group_indices()
    #' Get the median values for each clinical variable 
    #' Usually the variables are the top related to the factor used for clustering 



      # Get the median per cluster and scale it for the ehatmap 
      col_data_t<-tibble(col_data)


      means_by_cluster <- col_data_t %>% 
            dplyr::select(-PATNO) %>%
            group_by(cluster) %>%
            mutate_if(is.character, as.numeric) %>%
        summarise_all( funs(median(., na.rm = TRUE)))%>%
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

      means_by_cluster_na
      jpeg(outfile_cl_heatmap, width=7, height=3, units='in', res=200)

        p<-pheatmap(means_by_cluster_na, 
        labels_row= seq(1:k_centers_m),
        clustering_method = 'centroid',
        color = colorRampPalette(c("turquoise4", "white", "red"))(50))
        plot(p)
      dev.off()
      return(p)


}
plot_heatmap_median_by_cluster(col_data)

col_data<-samples_metadata(MOFAobjectPD_sel)[c(colnames(estimations),clust_name, 'PATNO')]
col_data$cluster<-col_data[, clust_name]; col_data[, clust_name]<-NULL
col_data$NK
outfile_cl_heatmap<-paste0(cluster_params_dir, '/heatmap_means_cells_new' ,  '.png')
outfile_cl_heatmap
plot_heatmap_median_by_cluster(col_data)



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


#install.packages('reshape')
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































