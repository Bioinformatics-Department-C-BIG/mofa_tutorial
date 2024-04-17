

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

DIFF_VAR= 'NP2PTOT_LOG'
DIFF_VAR= 'moca'


prot_de_mode<-'u'


tissue_un<-'Cerebrospinal Fluid'
tissue_un<-'Plasma'

prot_de_mode = 't'
TISSUE='Plasma';tissue_un<-'Plasma'
tissue ='Plasma'; tissue_un<-'Plasma'
DIFF_VAR= 'NP3TOT_LOG'


#source(paste0(script_dir, 'ppmi/extract_metadata.R'))
#source(paste0(script_dir, 'ppmi/analyse_clinical_vars.R'))


# Estimation of cell types from VSN rna data
# source(paste0(script_dir, 'ppmi/estimate_cell_types.R'))


#source(paste0(script_dir, 'ppmi/utils.R'))

source(paste0(script_dir, 'ppmi/mofa_application_ppmi_all_visits.R'))
source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))


# Plotting for mofa run not needed for clustering analysis, excluded  

#source(paste0(script_dir, 'ppmi/mofa_analysis_plots.R'))
#source(paste0(script_dir, 'ppmi/mofa_enrich.R')) # SET TO FALSE IF EXISTS? 

#source(paste0(script_dir, 'ppmi/mofa_enrich_pgcse.R')) # SET TO FALSE IF EXISTS? 


source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_analysis.R'))
# source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_plots.R'))

#source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' ))

 
formula_deseq_format='n'

cell_corr_deseq<-TRUE


VISIT_COMP='V08'
#VISIT='V08'
#prefix='mirnas_'
TISSUE='CSF'
TISSUE='Plasma'


#DIFF_VAR = 'NP2PTOT_LOG'


DIFF_VAR = 'NP2PTOT_LOG'
DIFF_VAR = 'moca'

diff_vars<-c( 'NP2PTOT_LOG', 'moca' ,'NP2PTOT_LOG')
diff_vars<-c( 'NP2PTOT_LOG','sft', 'NP3TOT_LOG',  'moca')

DIFF_VAR= 'updrs3_score_on_LOG'
DIFF_VAR = 'sft'



formula_deseq_format='n' # so far e only run all for mirnas 
formula_deseq_format='all' # so far we only run all for mirnas 

formula_deseq_format = 'age'
formula_deseq_format='age' # so far we only run all for mirnas 
formula_deseq_format='n' # so far we only run all for mirnas 

ONT='BP'
cell_corr_deseq = TRUE;
process_mirnas = FALSE; 
cell_corr_deseq=FALSE
VISIT_COMP='BL'

source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))


vis_comps<-c('V08',  'V06','V04', 'BL')
#vis_comps<-c('V08', 'BL')

VISIT_COMP='V04'
#vis_comps<-c('V08')

sig_only =FALSE
DIFF_VAR

run_all=TRUE 
if (run_all){
    for (DIFF_VAR in c(diff_vars)){
        print(DIFF_VAR)
        for (cell_corr_deseq in c(FALSE,  TRUE)){


                for (VISIT_COMP in vis_comps){
                print(VISIT_COMP)
                    source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))

                }


        }
        source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' ))


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
diff_vars<-c( 'NP2PTOT_LOG')

source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' ))



tissue_un<-'Plasma'




prot_de_mode<-'u'

#tissue_un<-'Plasma';tissue ='Plasma';

diff_vars<-c('NP3TOT_LOG')

DIFF_VAR= 'NP3TOT_LOG'
DIFF_VAR= 'NP2PTOT_LOG'

TISSUE='Plasma';

TISSUE='Plasma';


TISSUE='Plasma'
tissue_un<-'Plasma'
prot_de_mode = 't';

tissue_un<-'Cerebrospinal Fluid';tissue<-'Cerebrospinal Fluid'; 

visit_comps = c('V06', 'BL', 'V04', 'V08')
visit_comps  = c('V08' )
DIFF_VAR;tissue_un;prot_de_mode
sig_only =FALSE
VISIT_COMP='V06'
DIFF_VAR

# targeted plasma: not enough samples from cluster 1
    for (DIFF_VAR in c(diff_vars)){

        for (VISIT_COMP in visit_comps){
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

























































