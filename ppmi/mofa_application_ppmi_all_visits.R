source(paste0('ppmi/setup_os.R'))
#BiocManager::install('GOSemSim')
# SCENARIOS: 
# select cohort: 1,2,3,4: PD, Prodromal, , Healthy Control
# select visit: ALL, V02, V04, V06, V08 
library(MOFA2)
#install.packages("broom", type="binary")

library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)
library('MultiAssayExperiment')
source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))
source(paste0(script_dir,'ppmi/mofa_utils.R'))
source(paste0(script_dir,'ppmi/utils.R'))


split=FALSE
run_rna_mirna=FALSE
if (split){
  N_FACTORS=8
}
VISIT=c('V08');
run_vsn=TRUE
## tissue is set in the config
use_signif=FALSE
process_mirnas=FALSE
run_mofa_complete<-FALSE

source(paste0(script_dir, '/ppmi/config.R'))
source(paste0(script_dir, '/ppmi/mofa_config.R'))

# metadata source 
metadata_output<-paste0(output_files, 'combined.csv')
combined_all_original<-read.csv2(metadata_output)
metadata_output<-paste0(output_files, 'combined_log.csv')
combined_bl<-read.csv2(metadata_output)
scale_views=TRUE


### preferable settings based on the scenario 
## MORE samples, more factors ! 
# TODO: need a better way to decide how to run this 
if (run_mofa_complete){
  N_FACTORS=9
}else{
  #N_FACTORS=12
  N_FACTORS=15  ### so far gives the best of the corelations for TOP_GN=0.30, MN=0.5, PN=0.9
}






run_mofa_get_cors<-function(N_FACTORS){
  ### run mofa and write stats of corelation to file!! 
  #'
  #'

  ## histograms to check normality pattern 
  #create_hist(data_full['RNA'], 'RNA')
  #create_hist(data_full['miRNA'], 'miRNA')
  
  
  ## get list of three mats 
  data_full=prepare_multi_data(p_params, param_str_g, param_str_m, mofa_params)
  # create multi experiment 
  mofa_multi<-create_multi_experiment(data_full, combined_bl)
  mofa_multi_complete_all<-mofa_multi[,complete.cases(mofa_multi)]
  
  mofa_multi_rna_mir<-subsetByAssay(mofa_multi, c('RNA', 'miRNA'))
  mofa_multi_rna_mir_complete<-mofa_multi_rna_mir[,complete.cases(mofa_multi_rna_mir)]
  
  
  
  mofa_multi_to_use<- if(run_mofa_complete){mofa_multi_to_use=mofa_multi_complete}else{
    mofa_multi_to_use=mofa_multi
  }
  MOFAobject=create_mofa(mofa_multi_to_use)
  dir.create(outdir, showWarnings = FALSE)
  MOFAobject<-run_mofa_wrapper(MOFAobject, outdir )
  
  
  
  
  ##### Basic stats
  
  
  cors_both<-get_correlations(MOFAobject, c('CONCOHORT'))
  cors_pearson=cors_both[[2]]
  cors_t<-paste(round(cors_pearson[,'CONCOHORT'], digits=2), collapse=', ')
  max_cor<-round(max(cors_pearson), digits=2)
  print(cors_t)
  df_stats=  c( TOP_PN, TOP_GN, MIN_COUNT_G, TOP_MN, MIN_COUNT_M, mofa_params, sel_coh_s,VISIT_S,  scale_views[1],  use_signif,
                run_mofa_complete, N_FACTORS,cors_t , max_cor )
  
  write.table(t(df_stats), paste0(outdir_orig,'all_stats.csv'), append=TRUE,sep=',', col.names = FALSE)
  return(MOFAobject)
}


# n_factors best=15
for (N_FACTORS in c(15)){
  ## MOFA parameters, set directory 
  mofa_params<-paste0(N_FACTORS,'_sig_',  use_signif,'complete', run_mofa_complete )
  out_params<- paste0( 'p_', p_params, 'g_', g_params, 'm_', m_params, mofa_params, '_coh_', sel_coh_s,'_', VISIT_S, '_', scale_views[1])
  outdir = paste0(outdir_orig,out_params, '_split_', split , '/');outdir
  dir.create(outdir, showWarnings = FALSE)
  MOFAobject=run_mofa_get_cors(N_FACTORS)
}














