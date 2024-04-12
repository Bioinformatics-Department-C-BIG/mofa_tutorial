
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
clust_vars<-c('NP2PTOT_LOG', 'moca', 'NP3TOT_LOG', 'scopa', 'updrs3_score_on', 'sft')

get_factors_for_metric('updrs3_score_on')
facet_rows = 2





#### Boxplots ####
clust_vars=c('NP3TOT_LOG')

y_clust = 'NP3TOT_LOG'

#sapply(clust_vars, function(y_clust){

  for (y_clust in  clust_vars){
  clust_name = paste0(y_clust, '_clust')
  ## check if there are clusters for this variable





  #clust_name %in% colnames(met)
  if (clust_name %in% colnames(met)){
    # for each cluster create boxplot 


    fact = get_factors_for_metric(y_clust)


    # also write vars for each cluster 

    fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))

    cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors))
    cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );


          #### medians write csv ####
          # 1. anova -  which are the significant variables? 
        library(rstatix)
        diff_variables_box= c('NP3TOT','NP3TOT_V14','NP2PTOT', 'sft', 'moca', 'AGE', 'scopa', 'Neutrophil.Lymphocyte', 'SEX')
        col_data<-samples_metadata(MOFAobject_sel)[c(diff_variables_box,clust_name,'PATNO', colnames(estimations))]



        col_data$cluster<-col_data[, clust_name]; col_data[, clust_name]<-NULL
        col_data$cluster[col_data$cluster %in% c('HC')]<-0
        col_data$cluster<-as.numeric(col_data$cluster)


          col_data_t<-tibble(col_data)
          col_data_t_melt <- reshape2::melt(col_data_t, id.vars=c('cluster'))
          col_data_t_melt_not_na<-col_data_t_melt[!is.na(col_data_t_melt$value),]
        col_data_t_melt_not_na[is.na(col_data_t_melt_not_na) | col_data_t_melt_not_na=="Inf"] = NA
        col_data_t_melt_not_na<- na.omit(col_data_t_melt_not_na)




      ## ANALYSIS OF VARIANCE ### find out the sig varibles 
      variable = 'T.gd.Vd2'
      #' @param aov_sig_names the variables that are sig between the 4 groups  - including controls!! 
      #' pass on to the box plots etc. 

      aov_res_ll<-sapply(c(diff_variables_box, colnames(estimations)),function(variable){
        print(variable)
              col_data_t[c('cluster', variable)]
              col_data_t_1<-col_data_t[c('cluster', variable)]
              col_data_t_1[,'var'] = col_data_t_1[,variable]
              col_data_t_1

              tryCatch({
                na.omit(col_data_t_1[c('cluster', 'var')]$var)
                    res.aov<-col_data_t_1[c('cluster', 'var')] %>%
              anova_test(cluster ~ var) 
              res.aov
              }, 
              error = function(cond){ NA}
              )
          
              return(res.aov)
            }
)

# Choose boxplots based on AOV? 
aov_sig<-t(data.frame(aov_res_ll)[c('p', 'p<.05'),])
print(aov_sig)
aov_sig_names=rownames(aov_sig[aov_sig[,'p<.05'] == '*',])
print(aov_sig_names)



write.csv(t(aov_res_ll), paste0(cluster_params_dir,'/aov.csv'))


cluster_params_dir

## wilcox - compare with controls each cluster ##
# medians
medians_fname = paste0(cluster_params_dir, '/medians.csv')
medians_fname_tex = paste0(cluster_params_dir, '/medians.tex')

medians_cells_fname= paste0(cluster_params_dir, '/medians_cells.csv')

DataControl  = col_data %>%
            filter(cluster==0) 
DataControl

 var_x = 'moca'
pvals_all<-lapply(c(diff_variables_box, colnames(estimations)), function(var_x){
  DataControl$Value = DataControl[, var_x]
  col_data$Value = col_data[,var_x]
  col_data<-as.data.frame(apply(col_data,2,as.numeric))
  col_data$cluster
  DataControl$moca

  WILC<-col_data %>% tibble() %>%
  filter(cluster!=0) %>%
  group_by(cluster)
  
  
  pvals<-tryCatch({WILC %>% summarise(p_value = wilcox.test(DataControl$Value, Value, exact = FALSE)$p.value)}, 
  error = function(cond){
    NA
  })
  return(pvals)


}
)
names(pvals_all)<-c(diff_variables_box, colnames(estimations))
pvals_all_1<-pvals_all[!is.na(pvals_all)]
pvals_all_1
pvals_sig_any<-lapply(pvals_all_1, function(x){

    #print(as.data.frame(x[[2]]))
    print(any(x[[2]]<0.05))
    #print(x[[2]])

  
  
})
pvals_sig_any_true<-pvals_sig_any[unlist(pvals_sig_any)& !is.na(unlist(pvals_sig_any))]
pvals_sig_any_true_names<-names(pvals_sig_any_true)


          means_by_cluster <- col_data_t %>% 
          #      dplyr::select(-PATNO) %>%
                group_by(cluster) %>%
                mutate_if(is.character, as.numeric) %>%
            summarise_all( funs(median(., na.rm = TRUE)))%>% t()%>%
            as.data.frame() # %>%  dplyr::select(-cluster) %>%
          # as.data.frame()

    #      means_by_cluster[diff_variables_box,]<-round(means_by_cluster[diff_variables_box,])

          
          sex_means<-t(col_data_t %>% 
          group_by(cluster) %>% 
            mutate_if(is.character, as.numeric) %>%
                summarize_all(mean)) %>%  t()%>%
            as.data.frame()*100 
          

         write.csv(round(means_by_cluster[aov_sig_names,]),medians_fname )
         write.csv(means_by_cluster[colnames(estimations),],medians_cells_fname )
         rbind(means_by_cluster[aov_sig_names,], SEX= round(sex_means$SEX ))
         means_by_cluster['SEX',]= round(sex_means$SEX )

         round(means_by_cluster[diff_variables_box,], digits=0)
          print(xtable( format(means_by_cluster[c(aov_sig_names,'AGE','SEX'),], digits=1),
          type = "latex", 
          file = medians_fname_tex ))

          # separated to print with 2 digits
        aov_sig_names_cells<-aov_sig_names[aov_sig_names %in% colnames(estimations)]
            print(xtable( format(means_by_cluster[c('Neutrophil.Lymphocyte',aov_sig_names_cells),], digits=2),
          type = "latex", 
          file = medians_fname_tex ))






# 1. medians
# 2. get signifince  with zero cluster - wilcoxon
# 3. 
          IQR_by_cluster <- col_data_t %>% 
                dplyr::select(-PATNO) %>%
                group_by(cluster) %>%
                mutate_if(is.character, as.numeric) %>%
            summarise_all( funs(IQR(., na.rm = TRUE)))%>%
            as.data.frame() # %>%  dplyr::select(-cluster) %>

    }else{
      print(paste0(clust_name, '_clust does not exist'))
    }



  #### WRITE BOXPLOTS FOR THE SIG ONLY ####

  
    variates_to_p<-get_covariates_cells(y_clust, thresh=8)
    variates_to_p = c('B.Naive', 'B.Memory', 'Basophils.LD', 'Monocytes.C', 'Neutrophils.LD', 'T.CD4.Memory', 'T.CD4.Naive', 'T.CD8.Memory', 
    'T.CD8.Naive', 'Neutrophils....'  )
    variates_to_p<-variates_to_p[variates_to_p%in% colnames(met)]



    fact=get_factors_for_metric(y_clust)


    # also write vars for each cluster 
    write_vars_output(MOFAobject, vars_by_factor, factors=fact)

    fact_s=paste(fact[order(fact)], collapse='_'); print(paste(y_clust, fact_s))

    cluster_params<-paste0(fact_s ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors))
    cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );

    bn_all_fname<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,  '.png')

    dir.create(cluster_params_dir, recursive=TRUE)
    diff_variables_to_p<-diff_variables_to_p[diff_variables_to_p %in% colnames(met)] 


    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,  c(diff_variables_to_p), width=8+length(c(diff_variables_to_p))/facet_rows, 
    height=1+1.5*facet_rows, bn=bn_all_fname, facet_rows = 1, 
    text='')

    file.create(paste0(cluster_params_dir, '/', clust_name ))

    bn_all_fname<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,'cell_types',  '.png')
    variates_to_p = colnames(estimations)[!colnames(estimations) %in% c('T.gd.Vd2',"Plasmablasts" , 'mDCs', 'Monocytes.NC.I' )] # TODO: check for zero variance to exclude 
    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,  c(variates_to_p), width=8+length(c(variates_to_p))/facet_rows , bn=bn_all_fname, facet_rows = facet_rows, 
    text='')


    variates_to_p = intersect(pvals_sig_any_true_names, aov_sig_names_cells)

    bn_all_fname_sig<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,'_sig',  '.png')
    variates_to_p  = aov_sig_names[!aov_sig_names %in% c(colnames(estimations), 'scopa', 'NP3TOT_V14')]
    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,  c(variates_to_p), width=4+length(c(variates_to_p))/facet_rows ,
     bn=bn_all_fname_sig, facet_rows = facet_rows,     height=4,

    text='')

    variates_to_p =  aov_sig_names_cells # TODO: check for zero variance to exclude
    variates_to_p = intersect(pvals_sig_any_true_names, aov_sig_names_cells)
  variates_to_p = variates_to_p[!variates_to_p %in% c('Basophils.LD', 'T.CD8.Naive', 'T.gd.non.Vd2')]
  facet_rows = 1
    bn_all_fname_sig<-paste0(cluster_params_dir, '/all_vars_g_' ,sel_group,'_cells_sig',  '.png')
    boxplot_by_cluster_multiple(met=met, clust_name=clust_name,  c(variates_to_p), width=4+length(c(variates_to_p))/facet_rows ,
    height=4,
     bn=bn_all_fname_sig, facet_rows = facet_rows, 
    text='', plot_box=FALSE, add_caption = FALSE)






  }



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

dir.create(cluster_params_dir, recursive=TRUE)
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







plot_heatmap_median_by_cluster<-function(col_data){
      #
    #means_by_cluster %>% group_indices()
    #' Get the median values for each clinical variable 
    #' Usually the variables are the top related to the factor used for clustering 



      # Get the median per cluster and scale it for the heatmap 
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




## Split Mofa objects by cluster 
diff_var='sft'
clust_name<-paste0(diff_var,'_clust')
clust_id=3
fact<-get_factors_for_metric(diff_var)

graphics.off()
for (clust_id in c(1,2,3)){


      sm<-MOFAobject_sel@samples_metadata
      sel_group<-sm$PATNO_EVENT_ID[sm[, clust_name] %in% c(clust_id)]
      sel_group

      MOFAobject_cl<-subset_samples(MOFAobject_sel, sel_group)

      plot_data_overview(MOFAobject_cl)
      do_name<-paste0(outdir, '/clustering/',paste0(fact_s, collapse='_'),  '/data_overview' ,clust_id, '.jpeg' )
      jpeg(do_name)
      pl<-plot_data_overview(MOFAobject_cl)
      show(pl)
      dev.off()
}

do_name
































