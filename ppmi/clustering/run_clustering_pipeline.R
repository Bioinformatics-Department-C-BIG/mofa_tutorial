
script_dir<-paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/../../')
script_dir
cell_corr_deseq = TRUE
#VISIT='V08'
source(paste0('ppmi/setup_os.R'))
#source(paste0(script_dir, 'ppmi/setup_os.R'))
script_dir
source(paste0(script_dir, 'ppmi/mofa_application_ppmi_all_visits.R'))
source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))
#source(paste0(script_dir, 'ppmi/mofa_enrich.R')) # SET TO FALSE IF EXISTS? 



source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_analysis.R'))

VISIT_COMP='BL'
for (VISIT_COMP in c('BL', 'V06', 'V04', 'V08')){
    print(VISIT_COMP)
    source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))

}


# TODO: add mirs size effects using normalized data 
# correct for RIN and mapped bases 
## Add future scales 



