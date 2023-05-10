#install.packages('rstudioapi')

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))
print(script_dir)
suppressWarnings(library(sys))
suppressWarnings(library(data.table))


MIN_COUNT_G=100
MIN_COUNT_M=10


TOP_GN=0.10
TOP_MN=0.5
TOP_MN=0.75

### THE ONES THAT MOFA WILL USE! because to print we have set tgem in deseq_analysis
### if using signif


sel_coh=c(1,4)


sel_coh <- c(1,2)


sel_coh=c(1,2);
sel_subcoh=FALSE;
sel_subcoh=FALSE
sel_subcoh=c('INEXPD',  'INEXLRRK2', 'INEXSNCA');
sel_subcoh=c( 'INEXLRRK2', 'INEXSNCA');
sel_subcoh=c('INEXPD');

#1: INEXPD, INEXLRKK2, INEXSNCA 



VISIT_S=paste(VISIT,sep='_',collapse='-')
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')
if ( length(sel_subcoh)==1 && is.logical(sel_subcoh) && sel_subcoh==FALSE ){
  sel_subcoh_s=''
}else{
  sel_subcoh_s<-paste(sel_subcoh,sep='_',collapse='-')
}

sel_subcoh_s
g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 

param_str_m<-paste0('mirnas_',VISIT_S, '_', m_params ,'coh_',sel_coh_s, '_', sel_subcoh_s)
param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_', sel_subcoh_s)
param_str_m
param_str_m_f<-paste0('mirnas_', VISIT_S, '_',  MIN_COUNT_M, '_coh_',sel_coh_s, '_', sel_subcoh_s)
param_str_g_f<-paste0('rnas_', VISIT_S,  '_', MIN_COUNT_G,  '_coh_', sel_coh_s, '_' ,  sel_subcoh_s )
param_str_m

#### specific to rna seq 
output_1=paste0(data_dir, '/ppmi/output/')
output_files_orig<-output_1
output_files<-output_1
outdir_orig<-paste0(data_dir,'/ppmi/plots/')


### setup deseq formula 

formula_deseq<-'~AGE_AT_VISIT+SEX+EVENT_ID+COHORT'
formula_deseq2<-'~AGE_AT_VISIT+SEX+COHORT'
formula_deseq3<-'~PATNO+AGE_AT_VISIT+SEX'


des=gsub('~', '', formula_deseq2)


if (process_mirnas){
  prefix='mirnas_'
  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  param_str=param_str_m_f
  param_str_f=param_str_m_f
}else{
  prefix='rnas_'
   min.count=MIN_COUNT_G
   most_var=TOP_GN
   param_str=param_str_g_f
   param_str_f=param_str_g_f
}


input_file<-paste0(output_files, prefix, 'all_visits.csv')
# se file with all visits
se_file<-paste0(output_files, prefix,  'all_visits')
vsn_out_file<-paste0(output_files, param_str, '_vsn.csv')
deseq_file<-paste0(output_files, param_str_f, 'deseq.Rds')
outdir_s<-paste0(outdir_orig, '/single/', param_str_f, des)


######## PROCESS PROTEINS 
TISSUE='CSF'
TISSUE='Plasma'
TOP_PN=0.9


NA_PERCENT=0.9

NORMALIZED=TRUE;run_vsn=FALSE
run_vsn=TRUE
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')
VISIT_S=paste(VISIT,sep='_',collapse='-')

## VISIT_S to allow this to be more than one visits at once!! 
TOP_PN<-0.9

#### read in proteomics 
p_params_in<- paste0(  TISSUE, '_', NORMALIZED)
outdir_s_p<-paste0(outdir_orig, '/single/proteomics_', VISIT,'_norm_', substr(NORMALIZED,1,1),'vsn_', substr(run_vsn,1,1), '_coh_', sel_coh_s, '_',
                   sel_subcoh_s,  des, '/' )
p_params_FILE<- paste0(VISIT, '_', TISSUE, '_', NORMALIZED, '_',sel_coh_s,  sel_subcoh_s )
prot_vsn_se_filt_file<-paste0(output_files, p_params_FILE, '_vsn_se_filt.Rds')


## for mofa - run_vsn=TRUE





### mofa config 

#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#source(paste0(script_dir, '/config.R'))


