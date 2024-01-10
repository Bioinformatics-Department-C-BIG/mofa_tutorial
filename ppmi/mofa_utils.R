


get_correlations<-function(MOFAobject,covariates=c('CONCOHORT') ){
  cors<-correlate_factors_with_covariates(MOFAobject,
                                          covariates = covariates, 
                                          plot = "log_pval", 
                                          return_data = TRUE
                                          
  )
  
  cors_pearson<-correlate_factors_with_covariates(MOFAobject,
                                                  covariates = covariates, 
                                                  plot = "r", 
                                                  return_data = TRUE
                                                  
  )
  
  return(list(cors, cors_pearson))
}






correlate_factors_with_covariates_categ<-function(MOFAobject, covariates){
  #'
  #' performs Kruskal wallis test for mofa covariates 
  #' @param MOFAobject
  #' @param covariates the categorical covariates to test against
  #'
  #'
  NFACTORS<-MOFAobject@dimensions$K
   Z <- get_factors(MOFAobject, factors = 1:NFACTORS, 
                   as.data.frame = FALSE)[[1]]
  
  
  #Z <- do.call(rbind, Z)
  sm_pd<-MOFAobject@samples_metadata
  cat_covariates<-sm_pd[, covariates ]
  
  
  # Perform kruskal test for each factor and each covariate
  # Extract the pvalue 
  # TODO: consider adjusting? 
  cors_kruskal<-apply(Z, 2, function(z){
    
    apply(cat_covariates,2, 
          
          function(covariate){kruskal.test(z,
                                           covariate)$p.value })
    
  })
  
  
  
  
  cors_kruskal_log<-(-log10(cors_kruskal))
  th<-(-log10(0.05))
  #cors_kruskal_log[cors_kruskal_log < th ] <-0
  cors_kruskal_p<-cors_kruskal
  cors_kruskal_p[cors_kruskal_log < th ]<-0
  cors_kruskal_p
  #corrplot::corrplot(cors_kruskal_p)
  fname=paste0(outdir, '/covariates/cor_plot_categorical.jpeg'  )
  graphics.off()
  jpeg(fname, width=9, height=4, res=300, units='in')
  pheatmap(cors_kruskal_log, main='Correlations of factors with covariates \n Kruskal wallis test ')
        #    legend_labels=list(c('log10 p-adj.')))
  
  dev.off()
  return(cors_kruskal_log)
  
}




run_mofa_wrapper<-function(MOFAobject, outdir, force=FALSE, N_FACTORS=15 ){
  ### Run mofa and write to file
  #'
  #' @param MOFAobject 
  #' @param outdir
  #' 


  
  model_opts <- get_default_model_options(MOFAobject)
  data_opts <- get_default_data_options(MOFAobject)
  model_opts$num_factors <- N_FACTORS
  data_opts$scale_views=scale_views
  train_opts<-get_default_training_options(MOFAobject)

  
  
  MOFAobject <- prepare_mofa(MOFAobject,
                             model_options = model_opts,
                             data_options = data_opts, 
                             training_options = train_opts
  )
  
  mofa_file<-paste0(outdir,'mofa_ppmi.hdf5')
  
  if (file.exists(mofa_file) & !(force) ){
    pre_trained<-load_model(paste0(outdir,'mofa_ppmi.hdf5'))
    MOFAobject<-pre_trained
    
    
  }else {
    
    MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'), 
                           use_basilisk = TRUE)
  }
  
  
  return(MOFAobject)
  
}


#outdir1<-'D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/ppmi/plots/p_V08_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUE_split_FALSE - Copy/'
#pre_trained<-load_model(paste0(outdir,'mofa_ppmi.hdf5'))


factors=1

select_top_bottom_perc<-function(MOFAobject, view, factors, top_fr=.01 ){
  #'select top bottom features 
  #'#'
  #'factors
  #'factors
#  factors=1;view='RNA'
  ws<-get_weights(MOFAobject, views = view, factors=factors)[[1]]
  ws[order(ws)]
  cut_high<-top_fr; cut_low=1-top_fr
  high<-apply(ws,2, function(x){ ll<-as.data.frame(x) %>%
    top_frac(top_fr)
  return(rownames(ll))})
  low<-apply(ws,2, function(x){ ll<- -as.data.frame(x) %>%
    top_frac(top_fr)
  return(rownames(ll))})
  
  high_names<-reshape::melt(high)$value
  low_names<-reshape::melt(low)$value
  ws_union<-unique(c(high_names, low_names))
  return(ws_union)
}





concatenate_top_features<-function(MOFAobject, factors_all, view, top_fr){
  # get the top features from multiple factors and add the factor in second column
  f_all<-sapply(factors_all, function(f){
    select_top_bottom_perc(MOFAobject=MOFAobject,view=view, factors=f, top_fr=top_fr )
    }
    )
    colnames(f_all)<-factors_all
  f_features<-melt(f_all)[,2:3]
  colnames(f_features)<-c('Factor', 'feature')
  f_features$Factor
  return(f_features)
  
  
}





#object=MOFAobjectPD
#factors=c(6)
#cluster_samples_mofa_obj(MOFAobjectPD, k=2, factors=c(6))

cluster_samples_mofa_obj<-function(object, k, factors = "all", ...) 
{
       
        Z <- get_factors(object, factors = factors);
        if (is(Z, "list")) 
          Z <- do.call(rbind, Z)
        N <- nrow(Z)
        haveAllZ <- apply(Z, 1, function(x) all(!is.na(x)))
        if (!all(haveAllZ)) 
          warning(paste("Removing", sum(!haveAllZ), "samples with missing values on at least one factor"))
        Z <- Z[haveAllZ, ]
        
        Z_scaled <-apply(as.data.frame(Z), 2, scale )
        Z_scaled <-apply(as.data.frame(Z), 2, scale )
        
       # hist(Z_scaled[,1])
        #hist(Z_scaled[,2])
        
        kmeans.out <- kmeans(Z_scaled, centers = k, ...)
        return(kmeans.out)
      }


groups_from_mofa_factors<-function(patnos, MOFAobject_clusts, y_clust ){
  
  #' Obtain the molecular clusters from the mofa object 
  #' 
  #' @param MOFAobject description
  #' @param
  
  ### cluster by one factor 
  clust_name= paste0(y_clust, '_clust')
  sm_pd<-MOFAobject_clusts@samples_metadata;
  groups_kmeans<-sm_pd[, clust_name]
  pats<-sm_pd$PATNO
  kmeans_matched<-groups_kmeans[match(patnos, pats )]
  kmeans_grouping<-factor(kmeans_matched)
  kmeans_grouping
  return(kmeans_grouping)
  
}


#centers=3

###### CLUSTERS #### 


#rescale=TRUE
cluster_by_mofa_factors<-function(MOFAobject, factors,centers=2, rescale=FALSE, group=FALSE ){
  ###
  #' Cluster patients in a mofa object using the specified factors
  #' @param factors which factors to use for the clustering
  #' @return clusters_by_patno
  #' @
  #' 
  #' 
  gn<-length(MOFA2::groups_names(MOFAobject))

      all_groups_clusters<-sapply(1:gn, function(g) {
        Z <- get_factors(MOFAobject)[[g]]

        Z1 <- Z[,factors]
        if (rescale){
          Z1=scale(Z1)
        }
        groups_kmeans<-kmeans(Z1, centers=centers)
        return(groups_kmeans$cluster)
      }
      )
  
   groups_kmeans<-unlist( all_groups_clusters)


  return(groups_kmeans)
}




write_vars_output<-function(MOFAobject, vars_by_factor){
  #'
  #'
  #'write the variance plot 
  #' @param MOFAobject
  #' @param vars_by_factor
  #'
  #'
  vars_by_factor_f<-format(vars_by_factor*100, digits=2)
  write.table(format(vars_by_factor_f,digits = 2)
              ,paste0(outdir,'variance_explained.txt'), quote=FALSE)
  
  p3<-plot_variance_explained(MOFAobject, max_r2=20)+
    theme(axis.text.x=element_text(size=20), 
          axis.text.y=element_text(size=20))
  ggsave(paste0(outdir, 'variance_explained','.png'), plot=p3,
         width = 5, height=N_FACTORS/2,
         dpi=100)
}








boxplot_by_cluster<-function(met, clust_name, y, bn){
  #'
  #' Create boxplot by cluster 
  #' 
  #' @param met
  #' @param bn 
  #'  
  #' 
  clust_metric<-gsub('_clust', '', clust_name)
  
  met[,clust_metric ]<-as.numeric(met[,clust_metric])
  met<-met[!is.na(met[, clust_metric]),]
  print(paste('Using subset of  ', dim(met)[1], ' patients'))
  freqs<-paste0('n=', paste0(table(met[, clust_name]), collapse = ', '))
  
  
  #### PROPORTIONS OF BINARY VARS
  tot_med<-as.matrix(table(met[,c(clust_name, "PDMEDYN")])); paste_med<-paste0('Med: ' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ))
  tot_med<-as.matrix(table(met[,c(clust_name, "SEX")])); paste_sex<-paste('SEX:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  tot_med<-as.matrix(table(met[,c(clust_name, "PDSTATE")])); paste_state<-paste('PD state:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  tot_med<-as.matrix(table(met[,c(clust_name, "td_pigd")])); paste_tdpigd<-paste('td/pigd:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  tot_med<-as.matrix(table(met[,c(clust_name, "NHY")])); paste_nhy<-paste('H&Y:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  
  
  met[, clust_name]=as.factor(met[, clust_name])
  met[, y]=as.numeric(met[, y])
  
  
  k_centers<-max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)
  k_centers
  ## Add kruskal wallis to the total plot and separately to each one 
  ## 
  
  kw=NULL
  try(if (!all(is.na(met[, y]))){
    met_pd<-met[met$INEXPAGE%in% 'INEXPD',]
    kw<-kruskal.test(x=met_pd[, y], met_pd[, clust_name])
    
  })
  
  factors<-paste0(which(all_fs_diff[,clust_metric]), collapse=', ')
  
  my_comparisons=combn(levels(met[,clust_name]),2)
  a=max(met[, y], na.rm = TRUE)
  
  
  
  p<-ggplot(met ,aes_string(x=clust_name , y=y))+
    geom_boxplot(aes_string( x=clust_name, color=clust_name))+


    geom_signif(comparisons=list( c(1,2), c(2,3), c(1,3) ),
               aes_string(y=y), 
               y_position=c(a, a+0.5,a+1))+
            
    labs(title = paste(y),  
         subtitle=paste(freqs, '\n','Kruskal.wallis p.val', format(kw$p.value, digits=2)),
         caption = paste0('\n',
                          'factors: ',factors, '\n',
                          paste_med,  '\n',
                          paste_sex, '\n',
                          paste_state,  '\n',
                          paste_tdpigd, '\n', 
                          paste_nhy) )+
    theme(text =element_text(size=20))
  #plot.title = element_text(size = 30, color = "green")
    
  
  
  p
  print(bn)
  ggsave(bn, dpi=300)
  graphics.off()
  ## TODO: WILCOX TEST BY GROUP
  
  
}


#selected_covars=all_diff_variables_prog_conf
#MOFAobject_to_plot=MOFAobjectPD
#labels_col=FALSE

#selected_covars=all_diff_variables_prog_conf;
#factors = sel_factors_conf
#labels_col=TRUE
#height=1000
#res=200
# MOFAobject_to_plot=MOFAobjectPD_sel


get_labels<-function(selected_covars, labels_col=FALSE){
  # re-label selected clinical variables 
  if (labels_col){
    labels_col_plot<-mt_kv$V2[match(selected_covars,mt_kv$V1)]
    labels_col_plot[is.na(labels_col_plot)]<-selected_covars[is.na(labels_col_plot)]
  }else{
    labels_col_plot=selected_covars
  }
  return(labels_col_plot)
}


plot_covars_mofa<-function(selected_covars, fname, plot, factors,labels_col=FALSE, height=1000, 
                           MOFAobject_to_plot=MOFAobject, res=200){
  
  # filter if some do not exist in the colnames of metadata
  #apply(MOFAobjectPD@samples_metadata[,selected_covars3], 2, function(x) {length(which(duplicated(x)))==length(x)-1 })
  # first check if the requested names exist in the metadata 
  selected_covars=selected_covars[selected_covars %in% colnames(MOFAobject_to_plot@samples_metadata) ]
  
  
  sds<-apply(MOFAobject_to_plot@samples_metadata[,selected_covars], 2, sd, na.rm=TRUE)
  sd_na<-c(is.na(sds)|sds==0)
  
  # then check that the sd is not NA
  selected_covars=selected_covars[ !(sd_na) ]


  # remove if there are at least non-na values
  selected_covars = selected_covars[colSums(!is.na(MOFAobject_to_plot@samples_metadata[, selected_covars]))>3]
  
 P2_data<-correlate_factors_with_covariates(MOFAobject_to_plot,
                                        covariates =selected_covars , plot = plot,
                                        labels_col=get_labels(selected_covars, labels_col = TRUE), 
                                        factors = factors, 
                                        cluster_cols=TRUE, 
                                        return_data = TRUE)

  # keep only what is >0
  selected_covars = selected_covars[colSums(P2_data)>0]
  #plot='log_pval'

  
  jpeg(paste0(outdir,'/', fname,'.jpeg'), width = 1000+length(selected_covars)*20, height=height, res=res
         )

  # check which have zero corelation? 
  P2<-correlate_factors_with_covariates(MOFAobject_to_plot,
                                        covariates =selected_covars , 
                                        plot = plot,
                                        labels_col=get_labels(selected_covars, labels_col = TRUE), 
                                        factors = factors, 
                                        cluster_cols=TRUE)


  dev.off()
  
  
}
