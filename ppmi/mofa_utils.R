
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



run_mofa_wrapper<-function(MOFAobject, outdir, force=FALSE ){
  ### Run mofa and write to file
  #'
  #' @param MOFAobject 
  #' @param outdir
  #' 
  if (length(VISIT)>1){
    MOFAobject <- create_mofa(mofa_multi_complete, groups= mofa_multi_complete$EVENT_ID)
  }
  
  model_opts <- get_default_model_options(MOFAobject)
  data_opts <- get_default_data_options(MOFAobject)
  model_opts$num_factors <- N_FACTORS
  data_opts
  data_opts$scale_views=scale_views
  MOFAobject <- prepare_mofa(MOFAobject,
                             model_options = model_opts,
                             data_options = data_opts
  )
  
  mofa_file<-paste0(outdir,'mofa_ppmi.hdf5')
  if (file.exists(mofa_file) & !(force) ){
    pre_trained<-load_model(paste0(outdir,'mofa_ppmi.hdf5'))
    MOFAobject<-pre_trained
    
    
  }else {
    MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'), use_basilisk = TRUE)
  }
  
  return(MOFAobject)
  
}




