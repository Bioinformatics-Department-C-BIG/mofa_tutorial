


MIN_COUNT_G=100
MIN_COUNT_M=10
VISIT='BL'

VISIT=c('V04')
VISIT=('BL')
VISIT=('BL')




TOP_GN=0.1
TOP_MN=0.5


sel_coh=c(1,4)


sel_coh <- c(1,2)


sel_coh=c(1,2);





VISIT=c( 'V06')

VISIT_S=paste(VISIT,sep='_',collapse='-')
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')

g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 

param_str_m<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_')
param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  )



process_mirnas<-FALSE
if (process_mirnas){
  mirnas_file<-paste0(output_files, 'mirnas_all_visits.csv')
  mirnas_BL<-as.matrix(fread(mirnas_file, header=TRUE), rownames=1)
  
  raw_counts<-mirnas_BL
  
  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, param_str_m, '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')
  deseq_file<-paste0(output_files, param_str_m,'deseq.Rds')
  
  
}else{
  rnas_file<-paste0(output_files, 'rnas_all_visits.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
  
  raw_counts<-rnas_BL
  # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
  # raw_counts<-raw_counts %>% select(common_samples)
  
  min.count=MIN_COUNT_G
  most_var=TOP_GN
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, 'rnas_', param_str_g,  '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  
  highly_variable_outfile
  
  
  deseq_file<-paste0(output_files, param_str_g,'deseq.Rds')
  
  
}
deseq_file
raw_counts_all=raw_counts

