#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#source(paste0(script_dir, '/config.R'))
NORMALIZED=TRUE;
use_signif=FALSE
run_vsn=TRUE

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
  
  
  TOP_GN=0.10
  TOP_MN=0.3
  TOP_MN=0.25
  TOP_PN=0.75
  
  
}




### TODO move proteomics params here too!! 
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 
g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_') 

param_str_m<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_')
param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  )

### file specifically for mofa run 
p_params_out<- paste0(VISIT_S, '_',TISSUE, '_', TOP_PN, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)



highly_variable_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa.csv')
highly_variable_sign_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa_signif.csv')


