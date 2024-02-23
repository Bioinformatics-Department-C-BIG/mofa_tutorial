#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#source(paste0(script_dir, '/config.R'))
NORMALIZED=TRUE;
#use_signif=FALSE
run_vsn=TRUE
ruv=TRUE
scale_views=TRUE

split =FALSE

if (use_signif){
  TOP_GN=0.50
  TOP_MN=0.75
}else{
  TOP_GN=0.10
  TOP_MN=0.5
  TOP_MN=0.75

  
  ## MAKE A BIG ONE FOR Other purposes 
  TOP_MN=0.90
  TOP_GN=0.90
  

  ## MAKE A smaller one if we are using proteins too 
  TOP_MN=0.50
  TOP_GN=0.20
  
  
  
}
if (run_rna_mirna){
  ## MAKE A BIG ONE
  TOP_MN=0.90
  TOP_GN=0.30
  
}else{
  #  for (most_var in c(0.05, 0.1,0.2,0.3,  0.9,0.75,0.5)){
  
  ##
  TOP_GN=0.05
  TOP_MN=0.5
  
  TOP_GN=0.10
  TOP_MN=0.5
  

 
  

  
  ##mixomics
  TOP_MN=0.5
  TOP_GN=0.1
  TOP_PN=0.9
  
  if (run_mofa_complete){
    
    ### looks better with corelation to disease 
    
    TOP_GN=0.10
    TOP_MN=0.3
    TOP_MN=0.25
  }else{
    TOP_MN=0.50
    TOP_GN=0.30
    TOP_PN=0.9
    ## try with lower numbers for rnas...? 
    # But will GSEA work well?
    
    #TOP_MN=0.20
    #TOP_GN=0.10
    #TOP_PN=0.9
    
    
    ##mixomics
    TOP_MN=0.3
    TOP_GN=0.2
    TOP_PN=0.9
    
    ### OPTIMAL SO FAR 
    TOP_MN=0.50
    TOP_GN=0.30
    TOP_PN=0.9
    
    TOP_MN=0.50
    TOP_GN=0.30
    TOP_PN=0.9
    
    #
    TOP_MN=0.75
    TOP_GN=0.30
    TOP_PN=0.9
    
    #
    TOP_MN=0.50
    TOP_GN=0.40
    TOP_PN=0.9
    
    # FINAL
    TOP_MN=0.3
    
    TOP_MN=0.75
    TOP_MN=0.50
    
    TOP_GN=0.2
    TOP_GN=0.3
    TOP_GN=0.2
    
    
    TOP_PN=0.9
    
    
    if (use_signif){
      TOP_MN=0.75 # PLASMA
      TOP_MN=0.9# CSF PROTEOMICS
      TOP_GN=0.3
      TOP_PN=0.9
    }

  
  }

  
  
}

p_params_mofa<-p_params_plasma
