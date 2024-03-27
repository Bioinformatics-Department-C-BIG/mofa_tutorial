

#script_dir<-paste0(dirname(rstudioapi::getSourceEditorContext()$path), '/../')
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, 'ppmi/setup_os.R'))

# SCENARIOS: 
# select cohort: 1,2,3,4: PD, Prodromal, , Healthy Control
# select visit to run mofa: ALL, V02, V04, V06, V08 


library(MOFA2)
library(tidyverse)
library(psych)

library(data.table)
library(ggplot2)
library(ggpubr)
library(R.filesets)

library(dplyr)
library('MultiAssayExperiment')
source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))
source(paste0(script_dir,'ppmi/mofa_utils.R'))
source(paste0(script_dir, 'ppmi/predict_utils.R'))




split=FALSE
run_rna_mirna=FALSE
run_validation=FALSE
cell_corr_mofa=FALSE
cell_corr_deseq=TRUE
#if (split){
#  N_FACTORS=8
#}
VISIT=c('BL');




VISIT=c('BL','V08');
VISIT=c('BL','V04', 'V06',  'V08');
VISIT=c('V08');




run_vsn=TRUE
## tissue is set in the config
use_signif=FALSE
process_mirnas=FALSE
run_mofa_complete<-FALSE
run_rna_mirna<-FALSE
prot_mode = 'u'
#prot_mode = 't'

source(paste0(script_dir, '/ppmi/config.R'))
source(paste0(script_dir, '/ppmi/mofa_config.R'))
source(paste0(script_dir, '/ppmi/mofa_dirs.R'))
source(paste0(script_dir,'ppmi/utils.R'))

combined_bl_log<-load_metadata()
combined_bl <- combined_bl_log
decon<-loadRDS(paste0(output_files, '/decon.RDS'))
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

drop_factor_threshold = -1
drop_factor_threshold = -1
drop_factor_threshold = 0.001


print(TOP_PN)




run_mofa_get_cors<-function(mofa_multi_to_use, N_FACTORS, force=FALSE){
  ### run mofa and write stats of corelation to file!! 
  #'
  #'

 
  
  ##### Split and run many times ###
  ns<-dim(assays(mofa_multi_to_use)[['RNA']])[2]
  
  if (split){
    for (i in 1:ns){
      # separate and train on training data. Also save test data in outdir 
      mofa_multi_test<-mofa_multi_to_use[, ns]
      mofa_multi_train<-mofa_multi_to_use[, -ns]
      outdir=paste0(outdir, '_ns_', ns)
      
      if (complete.cases(mofa_multi_test)){
        # 
    
        MOFAobject_test=create_mofa(mofa_multi_test)
        
      }
      

    }
    }else{
      ## just use all data 
      mofa_multi_train=mofa_multi_to_use

    }
  
    MOFAobject=create_mofa(mofa_multi_train)
    if (length(VISIT)>1){
      # run mofa multi with event id as group
      MOFAobject <- create_mofa(mofa_multi_train, groups= mofa_multi_train$EVENT_ID)
      
    }
    
    
    dir.create(outdir, showWarnings = FALSE)
   # force=FALSE
    MOFAobject<-run_mofa_wrapper(MOFAobject, outdir, force=force, N_FACTORS=N_FACTORS, drop_factor_threshold=drop_factor_threshold )
    
    
  


 
  
  
  
  
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
#g_params

for (N_FACTORS in c(35)){
  ## MOFA parameters, set directory 
  #'  mofa_params<-paste0(N_FACTORS,'_sig_',  as.numeric(use_signif) ,'c_', as.numeric(run_mofa_complete)  )
  ruv_s<-(as.numeric(ruv))
   mofa_params<-paste0(N_FACTORS,'_sig_',  as.numeric(use_signif) ,'c_', as.numeric(run_mofa_complete)  )

  out_params<- paste0( 'p_',prot_mode, '_',TOP_PN,'_',p_params_plasma,'_', p_params_mofa, 'g_', g_params, 'm_', m_params, mofa_params, '_coh_', sel_coh_s,'_', VISIT_S, '_', 
                       as.numeric(scale_views[1]),'ruv_', as.numeric(ruv_s), '_c_',as.numeric(cell_corr_mofa), '_df_', drop_factor_threshold )
  
  outdir = paste0(outdir_orig,out_params, '_split_', as.numeric(split ));outdir
  dir.create(outdir, showWarnings = FALSE); outdir = paste0(outdir,'/' );outdir


#a_multi[,,4])

 ## get list of three matrices 
  # TODO: load all time points then filter? 

  data_full=prepare_multi_data(p_params, param_str_g_f, param_str_m_f,TOP_GN, TOP_MN,TOP_PN, mofa_params, prot_mode)
  # create multi-assay experiment object 
  mofa_multi<-create_multi_experiment(data_full, combined_bl_log)

  mofa_multi_complete_all<-mofa_multi[,complete.cases(mofa_multi)]


  # select if complete cases are used or samples with missing omics are included.
  # We need the multi assay experiment for other analysis too eg. time omics 
  mofa_multi_to_use<- if(run_mofa_complete){mofa_multi_to_use=mofa_multi_complete}else{
    mofa_multi_to_use=mofa_multi
  }

  MOFAobject=run_mofa_get_cors(mofa_multi_to_use, N_FACTORS, force=FALSE)
  
 
  
}


mofa_multi

inters_csf_rna<-intersectColumns(mofa_multi[, , c(2,3)])
inters_csf_rna$PATNO_EVENT_ID
mofa_multi[,,3]$PATNO_EVENT_ID %in% inters_csf_rna$PATNO_EVENT_ID
mofa_multi[,,3]$PATNO_EVENT_ID[!mofa_multi[,,3]$PATNO_EVENT_ID %in% inters_csf_rna$PATNO_EVENT_ID]


N_FACTORS =   MOFAobject@dimensions$K


length(intersect(colnames(assay(mofa_multi_to_use[,,2])),colnames(assay(mofa_multi_to_use[,,4])) ))
MOFAobject

MOFA2::plot_data_overview(MOFAobject)

## attach some extra clinical variables 
sel_sam=MOFA2::samples_names(MOFAobject)
length(sel_sam)
MOFA2::samples_metadata(MOFAobject)$group
meta_merged_ord=fetch_metadata_by_patient_visit(samples_metadata(MOFAobject)$PATNO_EVENT_ID)
length(meta_merged_ord$PATNO)
meta_merged_ord<-as.data.frame(meta_merged_ord)
meta_merged_ord$sample=meta_merged_ord$PATNO_EVENT_ID # MOFA needs a sample vector 
meta_merged_ord$group=MOFA2::samples_metadata(MOFAobject)$group # MOFA needs a sample vector 

#meta_merged_ord$sample=meta_merged_ord$PATNO_EVENT_ID
dim(samples_metadata(MOFAobject))
dim(as.data.frame(meta_merged_ord))
#MOFAobject@samples_metadata=as.data.frame(meta_merged_ord)
samples_metadata(MOFAobject)<-as.data.frame(meta_merged_ord)


sm<-samples_metadata(MOFAobject)


## ADD estimated cell types
estimations<-decon$proportions$qprogwc_ABIS_S0

colnames(estimations)

sm2<-cbind(samples_metadata(MOFAobject), estimations[match(sm$PATNO_EVENT_ID, rownames(estimations)),])


colnames(estimations) %in% colnames(samples_metadata(MOFAobject))

excl<-grepl( 'LAST_UPDATE|INFO_DT|TM|DT|ORIG_ENTRY|DATE|REC_ID|PAG_|SCENT_', colnames(sm2))
colnames(sm2)
sm3<-sm2[,!excl]
colnames(sm3)
samples_metadata(MOFAobject)<-as.data.frame(sm3)


## tests 
cors_real<-corr.test(sm2$Neutrophils.LD,sm2$Neutrophils....  )
cors_real$r

assay(mofa_multi_to_use[,,'proteomics_plasma'])
(as.data.frame(t(assay(mofa_multi_to_use[,,'proteomics_plasma']))))

rowVars(assay(mofa_multi_to_use[,,'proteomics_plasma']))
var(assay(mofa_multi_to_use[,,'proteomics_plasma'])['P07996',])

MOFAobject_orig<-MOFAobject

#mofa_multi_to_use[,, 'proteomics_csf']
































































































