#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#source(paste0(script_dir, '/config.R'))

use_signif=FALSE

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
  
  ## MAKE A BIG ONE
  TOP_MN=0.90
  TOP_GN=0.30
  
}


m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 
g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_') 

param_str_m<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_')
param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  )




highly_variable_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa.csv')
highly_variable_sign_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa_signif.csv')

