#install.packages('rstudioapi')

#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, 'ppmi/setup_os.R'))
print(script_dir)
suppressWarnings(library(sys))
suppressWarnings(library(data.table))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

MIN_COUNT_G=100
MIN_COUNT_M=10


TOP_GN=0.20
#TOP_MN=0.75
TOP_MN=0.50 # this is called inside mofa so check if it is the same as the mofa config... 

### THE ONES THAT MOFA WILL USE! because to print we have set tgem in deseq_analysis
### if using signif


# TODO: add in the config the ability to run BL with the V08 COMMON SAMPLES 
sel_coh=c(1,4)


sel_coh <- c(1,2)


sel_coh=c(1,2);

# ONLY PD 
sel_coh=c(4,2)


sel_coh=c(1);

sel_coh=c(1,2);
sel_coh=c(1,2);
sel_coh=c(1,2);


sel_subcoh=FALSE;
sel_subcoh=FALSE
sel_subcoh=c( 'INEXLRRK2', 'INEXSNCA');
sel_subcoh=FALSE
sel_subcoh=FALSE

sel_subcoh=c('INEXPD',  'INEXLRRK2', 'INEXSNCA');
sel_subcoh=c('INEXPD');

sel_subcoh=c('INEXPD',  'INEXLRRK2', 'INEXSNCA');
sel_subcoh=FALSE
sel_subcoh=c('INEXPD');

#cell_corr_deseq=TRUE

#sel_subcoh=c('INEXPD',  'INEXLRRK2', 'INEXSNCA');

#1: INEXPD, INEXLRKK2, INEXSNCA 




#VISIT=c( 'BL')



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
output_files<-paste0(data_dir, '/ppmi/output/')
outdir_orig<-paste0(data_dir,'/ppmi/plots/')


### setup deseq formula 

formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+COHORT'
formula_deseq2 = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+COHORT'
formula_deseq3<-'~PATNO+AGE_SCALED+SEX'
if (cell_corr_deseq) {
  formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+Neutrophil.Score+COHORT'
  
}else{
  formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+COHORT'
}

formula_deseq2 = formula_deseq


des=gsub('~', '', formula_deseq2)

process_mirnas


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



input_file<-paste0(output_files, prefix, 'all_visits.csv.gz')
input_file_mirs<-paste0(output_files, 'mirnas_', 'all_visits_norm.csv.gz')
input_file_mirs_norm<-input_file_mirs
# se file with all visits
se_file<-paste0(output_files, prefix,  'all_visits')
vsn_out_file<-paste0(output_files, param_str, '_vsn.csv')
vst_all_vis<-paste0(output_files, param_str, '_vst')
vst_cor_all_vis<-paste0(output_files, prefix, 'cell_corr_', cell_corr_mofa, 'vst_cor')  ## this is either filtered for rnas or full for mirnas  
vst_cor_all_vis


deseq_file<-paste0(output_files, param_str_f, 'deseq.Rds'); deseq_file
outdir_s<-paste0(outdir_orig, '/single/', param_str_f, des)


######## PROCESS PROTEINS #### load by tissue




TOP_PN_U=0.9




#TISSUE='Plasma'
#TISSUE='CSF'


NA_PERCENT=0.9

NORMALIZED=TRUE;run_vsn=FALSE
run_vsn=TRUE
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')
VISIT_S=paste(VISIT,sep='_',collapse='-')


## VISIT_S to allow this to be more than one visits at once!! 

TOP_PN<-0.5
TOP_PN
#### read in proteomics 
# this is untargeed
pr_project_id='9000'
p_params_in<- paste0( pr_project_id, '_',  TISSUE, '_', NORMALIZED)


outdir_s_p<-paste0(outdir_orig, '/single/proteomics_', VISIT,'_', pr_project_id, '_', TISSUE, '_norm_', substr(NORMALIZED,1,1),'vsn_', substr(run_vsn,1,1), '_coh_', sel_coh_s, '_',
                   sel_subcoh_s,  des, '/' )


p_params_FILE<- paste0(VISIT_S, '_', pr_project_id, '_',TISSUE, '_', NORMALIZED, '_',sel_coh_s,  sel_subcoh_s )



p_params<- paste0(VISIT_S, '_',pr_project_id, '_',TISSUE,  '_', substr(NORMALIZED,1,1), '_', sel_coh_s, sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)
p_params_csf<- paste0(VISIT_S, '_',pr_project_id, '_','CSF', '_', substr(NORMALIZED,1,1), '_', sel_coh_s, sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)
p_params_plasma<- paste0(VISIT_S, '_',pr_project_id, '_','Plasma', '_', substr(NORMALIZED,1,1), '_', sel_coh_s, sel_subcoh_s, 'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)


outdir_s_p<-paste0(outdir_orig, '/single/proteomics_', VISIT_S,'_',pr_project_id, '_',TISSUE, '_norm_', substr(NORMALIZED,1,1),'vsn_', substr(run_vsn,1,1), '_coh_', sel_coh_s, '_',
                   sel_subcoh_s,  des, '/' )


p_params_FILE<- paste0(VISIT_S, '_', pr_project_id, '_',TISSUE, '_', NORMALIZED, '_',sel_coh_s,  sel_subcoh_s )


prot_vsn_se_filt_file<-paste0(output_files, p_params_FILE, '_vsn_se_filt.Rds')


###### Untargeted proteins
 pr_un_project_id='177'


prot_untargeted_csf_f<-paste0(output_files, 'proteomics_',pr_un_project_id,'_','Cerebrospinal Fluid',  '.csv')
prot_untargeted_csf_vsn_f<-paste0(output_files, 'proteomics_',pr_un_project_id,'_','Cerebrospinal Fluid',  'vsn.csv')


prot_untargeted_plasma_f<-paste0(output_files, 'proteomics_',pr_un_project_id,'_','Plasma',  '.csv')
prot_untargeted_plasma_vsn_f<-paste0(output_files, 'proteomics_',pr_un_project_id,'_','Plasma', 'vsn.csv')



#tissue_un<-'Cerebrospinal Fluid'
#tissue_un<-'Plasma'

prot_untargeted_un_f<-paste0(output_files, 'proteomics_',pr_un_project_id,'_',tissue_un,  '.csv')
prot_untargeted_un_vsn_f<-paste0(output_files, 'proteomics_',pr_un_project_id,'_',tissue_un,  'vsn.csv')





























