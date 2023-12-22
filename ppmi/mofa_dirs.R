
### TODO move proteomics params here too!! 
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 
g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_') 

param_str_m_f<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_', sel_subcoh_s)
param_str_g_f<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  , sel_subcoh_s)

### file specifically for mofa run 
p_params_out<- paste0(VISIT_S, '_',TISSUE, '_', TOP_PN, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,sel_subcoh_s, 'vsn_', substr(run_vsn,1,1),
     'NA_', NA_PERCENT)


# TODO: DEFINE the deseq files 
#deseq_file_mirs<-
#  deseq_file_genes<-
#highly_variable_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa.csv')
#highly_variable_sign_outfile<-paste0(output_files, param_str,'_highly_variable_genes_mofa_signif.csv')


outdir_orig=paste0(data_dir,'ppmi/plots/')
output_files<- paste0(data_dir,'ppmi/output/')


