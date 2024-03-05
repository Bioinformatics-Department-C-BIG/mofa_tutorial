

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
# source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_plots.R'))



VISIT_COMP='BL'
VISIT_COMP='V08'

formula_deseq_format='n'

cell_corr_deseq<-TRUE

VISIT_COMP='BL'
VISIT_COMP='V06'
VISIT_COMP='BL'
VISIT_COMP='V08'
#VISIT='V08'
#prefix='mirnas_'
TISSUE='CSF'
TISSUE='Plasma'
DIFF_VAR = 'moca'
DIFF_VAR = 'NP2PTOT_LOG'

process_mirnas=TRUE 
process_mirnas = FALSE 

formula_deseq_format='n' # so far e only run all for mirnas 
formula_deseq_format='all' # so far we only run all for mirnas 
cell_corr_deseq=FALSE;

ONT='BP'
source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))

run_all=TRUE 
if (run_all){

for (cell_corr_deseq in c( FALSE)){


        for (VISIT_COMP in c('BL', 'V06', 'V04', 'V08')){
           print(VISIT_COMP)
            source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))

        }


}}

# after running all time points run time comparisons 

source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons_pathways_all_time.R'))



# for proteins only 
# TODO: change the config to allow plasma or csf 
# TISSUE=Plasma, CSF
#prefix='prot_'
 for (VISIT_COMP in c('V06', 'BL', 'V04', 'V08')){
#for (VISIT_COMP in c( 'V04', 'V08')){
    # compares each time point separatelty 
    source(paste0(script_dir, 'ppmi/clustering/DE_tutorial_ppmi_cluster_compare.R'))
    # ensure that the settings were not updated elsewhere
    print(paste0(TISSUE,' VISIT: ', VISIT_COMP ))
        }
# TODO: add mirs size effects using normalized data 
# correct for RIN and mapped bases 
## Add future scales 

# Run the proteins too 
# Concatenates all time points 
source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons_proteins_time.R'))


















