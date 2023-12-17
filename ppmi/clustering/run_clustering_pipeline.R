
script_dir<-paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/../')
script_dir
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, 'ppmi/setup_os.R'))

source(paste0(script_dir, 'ppmi/mofa_application_ppmi_all_visits.R'))
source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))
source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_analysis.R'))
source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))


source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))
