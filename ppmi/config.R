#install.packages('rstudioapi')

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))
print(script_dir)
library(sys)
library(data.table)


MIN_COUNT_G=100
MIN_COUNT_M=10



### if using signif
if (use_signif){
  TOP_GN=0.50
  TOP_MN=0.75
}else{
  TOP_GN=0.10
  TOP_MN=0.5
}


sel_coh=c(1,4)


sel_coh <- c(1,2)


sel_coh=c(1,2);


VISIT_S=paste(VISIT,sep='_',collapse='-')
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')

g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 

param_str_m<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_')
param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  )

param_str_m_f<-paste0('mirnas_',VISIT_S, '_',MIN_COUNT_M, '_coh_',sel_coh_s, '_')
param_str_g_f<-paste0('rnas_', VISIT_S, '_',MIN_COUNT_G, '_coh_', sel_coh_s, '_'  )


#### specific to rna seq 
output_1=paste0(data_dir, 'ppmi/output/')
output_files_orig<-output_1
output_files<-output_1
outdir_orig<-paste0(data_dir,'ppmi/plots/')


### setup deseq formula 

formula_deseq<-'~AGE_AT_VISIT+SEX+COHORT+EVENT_ID'
formula_deseq2<-'~AGE_AT_VISIT+SEX+COHORT'
formula_deseq3<-'~PATNO+AGE_AT_VISIT+SEX'


des=gsub('~', '', formula_deseq2)


des
process_mirnas<-TRUE
if (process_mirnas){
  input_file<-paste0(output_files, 'mirnas_all_visits.csv')
 
  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, param_str_m, '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_m_f,'_highly_variable_genes_mofa.csv')
  highly_variable_sign_outfile<-paste0(output_files, param_str_m_f,'_highly_variable_genes_mofa_signif.csv')
  
  deseq_file<-paste0(output_files, param_str_m, 'deseq.Rds')
  outdir_s<-paste0(outdir_orig, '/single/', param_str_m, des)
  
  
}else{
  input_file<-paste0(output_files, 'rnas_all_visits.csv')
  
  # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
  # raw_counts<-raw_counts %>% select(common_samples)
  
  min.count=MIN_COUNT_G
  most_var=TOP_GN
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, 'rnas_', param_str_g,  '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  highly_variable_sign_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa_signif', '.csv')
  
  highly_variable_outfile
  
  
  deseq_file<-paste0(output_files, param_str_g_f,'deseq.Rds')
  
  outdir_s<-paste0(outdir_orig, '/single/', param_str_g_f, des)
  
  
}


