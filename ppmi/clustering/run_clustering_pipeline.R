

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

#source(paste0(script_dir, 'ppmi/mofa_enrich_pgcse.R')) # SET TO FALSE IF EXISTS? 


source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_analysis.R'))
#source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_plots.R'))

#source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' ))



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


#DIFF_VAR = 'NP2PTOT_LOG'


DIFF_VAR = 'NP2PTOT_LOG'
DIFF_VAR = 'moca'

diff_vars<-c( 'NP2PTOT_LOG', 'moca','updrs3_score_on_LOG' ,'NP3TOT_LOG')
DIFF_VAR= 'updrs3_score_on_LOG'
DIFF_VAR = 'moca'
process_mirnas = TRUE 
process_mirnas = FALSE 


formula_deseq_format='n' # so far e only run all for mirnas 
formula_deseq_format='all' # so far we only run all for mirnas 

cell_corr_deseq=TRUE
formula_deseq_format='n' # so far we only run all for mirnas 

ONT='BP'
source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))


VISIT_COMP
vis_comps<-c('V08',  'V06','V04', 'BL')
#vis_comps<-c('V08')
cell_corr_deseq = TRUE
DIFF_VAR
run_all=TRUE 
if (run_all){
    for (DIFF_VAR in c(diff_vars)){
        print(DIFF_VAR)
        for (cell_corr_deseq in c( TRUE)){


                for (VISIT_COMP in vis_comps){
                print(VISIT_COMP)
                    source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))

                }


        }
        source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons_pathways_all_time.R'))

}
}

# after running all time points run time comparisons 




# for proteins only s
# TODO: change the config to allow plasma or csf 
# TISSUE=Plasma, CSF
#prefix='prot_'
#DIFF_VAR
diff_vars<-c( 'NP3TOT_LOG', 'NP2PTOT_LOG', 'moca')
diff_vars<-c( 'NP2PTOT_LOG', 'moca','updrs3_score_on_LOG' ,'NP3TOT_LOG')

DIFF_VAR= 'NP3TOT_LOG'

    for (DIFF_VAR in c(diff_vars)){

        for (VISIT_COMP in c('V06', 'BL', 'V04', 'V08')){
        #for (VISIT_COMP in c( 'V04', 'V08')){
            # compares each time point separatelty 
            source(paste0(script_dir, 'ppmi/clustering/DE_tutorial_ppmi_cluster_compare.R'))
            # ensure that the settings were not updated elsewhere
            print(paste0(tissue,' VISIT: ', VISIT_COMP ))
            print(fact)
            print(cluster_params_dir)
                }
    source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons_proteins_time.R'))

    }
# TODO: add mirs size effects using normalized data 
# correct for RIN and mapped bases 
## Add future scales 

# Run the proteins too 
# Concatenates all time points 




































