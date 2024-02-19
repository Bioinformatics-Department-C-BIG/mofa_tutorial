

cell_corr_deseq = TRUE
#VISIT='V08'
print(paste('wd: ',getwd()))

try(
    source(paste0('ppmi/setup_os.R'))
)

try(
   source(paste0('../../ppmi/setup_os.R')), silent = TRUE

)

#source(paste0(script_dir, 'ppmi/setup_os.R'))

source(paste0(script_dir, 'ppmi/mofa_application_ppmi_all_visits.R'))
source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))


# Plotting for mofa run not needed for clustering analysis, excluded  
#source(paste0(script_dir, 'ppmi/mofa_analysis_plots.R'))
#source(paste0(script_dir, 'ppmi/mofa_enrich.R')) # SET TO FALSE IF EXISTS? 



source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_analysis.R'))
source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_plots.R'))

VISIT_COMP='BL'
VISIT_COMP='V08'

formula_deseq_format='n'

cell_corr_deseq<-TRUE
cell_corr_deseq<-FALSE

VISIT_COMP='BL'
VISIT_COMP='V06'
VISIT_COMP='BL'
VISIT_COMP='V08'
formula_deseq_format='all'

for (cell_corr_deseq in c(TRUE, FALSE)){


        for (VISIT_COMP in c('BL', 'V06', 'V04', 'V08')){
        #   print(VISIT_COMP)
            source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))

        }
}

# TODO: add mirs size effects using normalized data 
# correct for RIN and mapped bases 
## Add future scales 




