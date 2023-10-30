


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




select_top_bottom_perc<-function(view, factors, top_fr=.01 ){
  #'select top bottom features 
  #'#'
  #'factors
  #'factors
#  factors=1;view='RNA'
  ws<-get_weights(MOFAobject, views = view, factors=factors)[[1]]
  print(factors)
  cut_high<-top_fr; cut_low=1-top_fr
  high<-apply(ws,2, function(x){ ll<-as.data.frame(x) %>%
    top_frac(top_fr)
  return(rownames(ll))})
  low<-apply(ws,2, function(x){ ll<--as.data.frame(x) %>%
    top_frac(top_fr)
  return(rownames(ll))})
  
  high_names<-melt(high)$value
  low_names<-melt(low)$value
  ws_union<-unique(c(high_names, low_names))
  return(ws_union)
}

#factors=heatmap_factors
concatenate_top_features<-function(factors, view, top_fr){
  
  f_all<-sapply(factors,function(f){
    select_top_bottom_perc(view=view, factors=f, top_fr=top_fr )
    }
    )
  f_features<-melt(f_all)[,2:3]
  colnames(f_features)<-c('Factor', 'feature')
  
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







