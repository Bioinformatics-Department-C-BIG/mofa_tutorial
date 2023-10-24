

#script_dir<-paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/../')
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, 'ppmi/setup_os.R'))

#source(paste0('/Users/efiathieniti/Documents/GitHub/mofa_tutorial/ppmi/setup_os.R'))
#BiocManager::install('GOSemSim')
# SCENARIOS: 
# select cohort: 1,2,3,4: PD, Prodromal, , Healthy Control
# select visit: ALL, V02, V04, V06, V08 
# BiocManager::install("MOFA2", version="1.8")
#etach('package:MOFA2', unload=TRUE)
#source("https://bioconductor.org/biocLite.R")


library(MOFA2)


library(data.table)
library(ggplot2)
library(ggpubr)
library(dplyr)
library('MultiAssayExperiment')
source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))
source(paste0(script_dir,'ppmi/mofa_utils.R'))
source(paste0(script_dir,'ppmi/utils.R'))
source(paste0(script_dir, 'ppmi/predict_utils.R'))


split=FALSE
run_rna_mirna=FALSE
run_validation=FALSE
#if (split){
#  N_FACTORS=8
#}
VISIT=c('BL');
VISIT=c('BL','V04', 'V06',  'V08');
VISIT=c('BL','V08');
VISIT=c('V08');


run_vsn=TRUE
## tissue is set in the config
use_signif=FALSE
process_mirnas=FALSE
run_mofa_complete<-FALSE

source(paste0(script_dir, '/ppmi/config.R'))
source(paste0(script_dir, '/ppmi/mofa_config.R'))
source(paste0(script_dir, '/ppmi/mofa_dirs.R'))

# metadata source 
metadata_output<-paste0(output_files, 'combined.csv')
combined_all_original<-read.csv2(metadata_output)
metadata_output<-paste0(output_files, 'combined_log.csv') 
combined_bl_log<-read.csv2(metadata_output) # combined_bl_log holds the updated data , log, scaled, future visits 
combined_bl_log$GBA_PATHVAR

combined_bl<-combined_all_original
combined_bl<-combined_bl_log

combined_bl_log$RBD_TOT
all(is.na(combined_bl_log$updrs2_score_BL))

#combined_all_original[combined_all_original$PATNO_EVENT_ID %in% sel_sam, ]$CONCOHORT
#combined_bl_log[combined_bl_log$PATNO_EVENT_ID %in% sel_sam, ]$CONCOHORT

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









run_mofa_get_cors<-function(N_FACTORS, force=FALSE){
  ### run mofa and write stats of corelation to file!! 
  #'
  #'

  ## histograms to check normality pattern 
  #create_hist(data_full['RNA'], 'RNA')
  #create_hist(data_full['miRNA'], 'miRNA')
  
  #combined_bl=combined_bl_log
  ## get list of three mats 
  data_full=prepare_multi_data(p_params, param_str_g_f, param_str_m_f,TOP_GN, TOP_MN, mofa_params)
  # create multi experiment 
  mofa_multi<-create_multi_experiment(data_full, combined_bl)
  mofa_multi_complete_all<-mofa_multi[,complete.cases(mofa_multi)]
  
  #mofa_multi_rna_mir<-subsetByAssay(mofa_multi, c('RNA', 'miRNA'))
 # mofa_multi_rna_mir_complete<-mofa_multi_rna_mir[,complete.cases(mofa_multi_rna_mir)]
  
  
  
  mofa_multi_to_use<- if(run_mofa_complete){mofa_multi_to_use=mofa_multi_complete}else{
    mofa_multi_to_use=mofa_multi
  }
  
  ##### Split and run many times ###
  ns<-dim(assays(mofa_multi_to_use)[['RNA']])[2]
  
  if (split){
    for (i in 1:ns){
      # separate and train on training data. Also save test data in outdir 
      mofa_multi_test<-mofa_multi_to_use[, ns]
      mofa_multi_train<-mofa_multi_to_use[, -ns]
      outdir=paste0(outdir, '_ns_', ns)
      
      if (complete.cases(mofa_multi_test)){
        
    
        MOFAobject_test=create_mofa(mofa_multi_test)
        
        
        
      }
      

    }
    }else{
      ## just use all data 
      mofa_multi_train=mofa_multi_to_use
      
    
    
    
    }
  
    MOFAobject=create_mofa(mofa_multi_train)
    #if (length(VISIT)>1){
    #  MOFAobject <- create_mofa(mofa_multi_train, groups= mofa_multi_train$EVENT_ID)
    ##  
    #}
    
    
    dir.create(outdir, showWarnings = FALSE)
   # force=FALSE
    MOFAobject<-run_mofa_wrapper(MOFAobject, outdir, force=force, N_FACTORS=N_FACTORS )
    
    
  


 
  
  
  
  
  ##### Basic stats
  #'
  #'
  #'
  
  
  if (length(sel_coh)>1){
    #' Only check correlations with cohort for more than one cohorts 
        cors_both<-get_correlations(MOFAobject, c('CONCOHORT'))
        cors_pearson=cors_both[[2]]
        cors_t<-paste(round(cors_pearson[,'CONCOHORT'], digits=2), collapse=', ')
        max_cor<-round(max(cors_pearson), digits=2)
        print(cors_t)
       
  
        if (length(VISIT)==1){
          
    
        ### write cross val
        ### todo: save random forest results 
        cors_pvalue=cors_both[[1]]
        sel_factors<-which(cors_pvalue>-log10(0.05))
        
        print(sel_factors)
        N_final<-MOFAobject@dimensions$K
        df_mofa <- as.data.frame(get_factors(MOFAobject, factors=1:N_final)[[1]])
        df_mofa$y<- as.factor(MOFAobject@samples_metadata$COHORT)
        df_mofa_age <- cbind(df_mofa,MOFAobject@samples_metadata[, c('AGE_SCALED', 'SEX')])
        
        if (run_validation){
          res_age_mofa<-run_train_validation( df=df_mofa_age)
          acc_mean<-mean(res_age_mofa[ 'Balanced Accuracy', na.rm=TRUE])
          print(acc_mean)
          
          df_stats=  c( TOP_PN, TOP_GN, MIN_COUNT_G, TOP_MN, MIN_COUNT_M, mofa_params, sel_coh_s,VISIT_S,  scale_views[1],  use_signif,
                        run_mofa_complete, N_FACTORS,cors_t , max_cor, acc_mean )
          
          write.table(t(df_stats), paste0(outdir_orig,'all_stats.csv'), append=TRUE,sep=',', col.names = FALSE)
        }
       
        }
    
  }
  
  
  return(MOFAobject)
  
  
}


# n_factors best=15
for (N_FACTORS in c(15)){
  ## MOFA parameters, set directory 
  #'
  mofa_params<-paste0(N_FACTORS,'_sig_',  use_signif,'complete', run_mofa_complete )
  out_params<- paste0( 'p_', p_params, 'g_', g_params, 'm_', m_params, mofa_params, '_coh_', sel_coh_s,'_', VISIT_S, '_', scale_views[1])
  outdir = paste0(outdir_orig,out_params, '_split_', split );outdir
  dir.create(outdir, showWarnings = FALSE)
  MOFAobject=run_mofa_get_cors(N_FACTORS, force=FALSE)
  outdir = paste0(outdir,'/' );outdir
  
  
  
}

## attach some extra clinical variables 
sel_sam=MOFAobject@samples_metadata$PATNO_EVENT_ID
length(sel_sam)
meta_merged_ord=fetch_metadata_by_patient_visit(MOFAobject@samples_metadata$PATNO_EVENT_ID)
length(meta_merged_ord$PATNO)
meta_merged_ord<-as.data.frame(meta_merged_ord)
meta_merged_ord$sample=meta_merged_ord$PATNO_EVENT_ID

MOFAobject@samples_metadata=as.data.frame(meta_merged_ord)
samples_metadata(MOFAobject)<-as.data.frame(meta_merged_ord)





