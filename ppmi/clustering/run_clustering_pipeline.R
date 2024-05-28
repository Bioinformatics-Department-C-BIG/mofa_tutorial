

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

# PREPROCESSING 

#source(paste0(script_dir, 'ppmi/proteomics/olink-preprocessing_ppmi_se.R'))
#source(paste0(script_dir, 'ppmi/proteomics/untargeted_preprocessing.R'))



source(paste0(script_dir, 'ppmi/mofa_application_ppmi_all_visits.R'))
source(paste0(script_dir, 'ppmi/mofa_analysis_time_diff.R'))


# Plotting for mofa run not needed for clustering analysis, excluded  

#source(paste0(script_dir, 'ppmi/mofa_analysis_plots.R'))
#source(paste0(script_dir, 'ppmi/mofa_enrich.R')) # SET TO FALSE IF EXISTS? 

#source(paste0(script_dir, 'ppmi/mofa_enrich_pgcse.R')) # SET TO FALSE IF EXISTS? 


source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_analysis.R'))
source(paste0(script_dir, 'ppmi/clustering/mofa_clustering_plots.R'))



formula_deseq_format='n'
cell_corr_deseq<-TRUE

#VISIT='V08'
#prefix='mirnas_'
TISSUE='CSF'
TISSUE='Plasma'


#DIFF_VAR = 'NP2PTOT_LOG'


DIFF_VAR = 'NP2PTOT_LOG'
DIFF_VAR = 'moca'

diff_vars<-c( 'NP2PTOT_LOG', 'moca' ,'NP3TOT_LOG')
diff_vars<-c('NP2PTOT_LOG','sft' ,'NP3TOT_LOG', 'moca' ,'NP2PTOT_LOG')
#diff_vars<-c('sft' )

DIFF_VAR= 'updrs3_score_on_LOG'
DIFF_VAR = 'moca'
DIFF_VAR = 'sft'
diff_vars<-c( 'COHORT')
diff_vars<-c( 'NP2PTOT_LOG','sft', 'NP3TOT_LOG',  'moca',  'NP3TOT_LOG_V14', 'NP3TOT_diff_V14', 'NP3TOT_V14', 'sft_V12' )
diff_vars<-c('NP3TOT_LOG' )
diff_vars<-c('NP3TOT_LOG')

DIFF_VAR = 'COHORT'
DIFF_VAR = 'NP3TOT_LOG'

formula_deseq_format='n' # so far e only run all for mirnas 
formula_deseq_format='all' # so far we only run all for mirnas 

cell_corr_deseq=FALSE
formula_deseq_format='n' # so far we only run all for mirnas 
formula_deseq_format='age' # so far we only run all for mirnas 

ONT='BP'
process_mirnas = FALSE 
source(paste0(script_dir, 'ppmi/config.R'))
VISIT_COMP='V08'

source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons.R'))
se_clust$kmeans_grouping



vis_comps<-c('V08',  'V06','V04', 'BL')
#vis_comps<-c('V08', 'BL')

VISIT_COMP='V08'
#vis_comps<-c('V08')
cell_corr_deseq = FALSE

sig_only =FALSE

 # NEED TO OBTAIN metrics from combined_bl_log_sel_mol
 source(paste0(script_dir,'/ppmi/clinical_variables_over_time.R' )) 
run_all=TRUE 
if (run_all){
    for (DIFF_VAR in c(diff_vars)){
        print(DIFF_VAR)
        for (cell_corr_deseq in c(  FALSE)){


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
diff_vars<-c( 'NP2PTOT_LOG', 'moca','NP3TOT_LOG')
diff_vars<-c( 'NP3TOT_LOG', 'NP2PTOT_LOG', 'moca', 'sft')
diff_vars<-c( 'NP3TOT_LOG')

DIFF_VAR= 'NP2PTOT_LOG'
DIFF_VAR= 'moca'




tissue_un<-'Plasma'




prot_de_mode<-'u'

#tissue_un<-'Plasma';tissue ='Plasma';

diff_vars<-c('sft', 'NP3TOT_LOG', 'NP2PTOT_LOG', 'moca')

DIFF_VAR= 'NP2PTOT_LOG'
DIFF_VAR= 'NP3TOT_LOG'

TISSUE='Plasma';

prot_de_mode = 'u';
TISSUE='Plasma';


TISSUE='CSF'
tissue_un<-'Plasma'
prot_de_mode = 'u';

tissue_un<-'Cerebrospinal Fluid';tissue<-'Cerebrospinal Fluid'; 
tissue_un<-'Plasma'
visit_comps  = c('V08' )
visit_comps = c('V06', 'BL', 'V04', 'V08')
visit_comps = c( 'V08','V06')

DIFF_VAR;tissue_un;prot_de_mode
sig_only =TRUE
VISIT_COMP='V08'

    for (DIFF_VAR in c(diff_vars)){

        for (tissue_un in c('Cerebrospinal Fluid','Plasma' )){
            for (VISIT_COMP in visit_comps){

           
        #for (VISIT_COMP in c( 'V04', 'V08')){
            # compares each time point separatelty 
            source(paste0(script_dir, 'ppmi/clustering/DE_tutorial_ppmi_cluster_compare.R'))
            # ensure that the settings were not updated elsewhere
            print(paste0(tissue,' VISIT: ', VISIT_COMP ))
            print(fact)
            print(cluster_params_dir)
                }

            source(paste0(script_dir, 'ppmi/clustering/cluster_comparisons_proteins_time.R')) # for each tissue combine times 

        }

    }
# TODO: add mirs size effects using normalized data 
# correct for RIN and mapped bases 
## Add future scales 

# Run the proteins too 
# Concatenates all time points 

















































































