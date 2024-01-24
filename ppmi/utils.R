
library(DESeq2)
library(edgeR)
library(org.Hs.eg.db)
library(sgof)
#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/factoextra")
library('factoextra')
library(plyr)
library(dplyr)
library('WGCNA') # needƒ empiricalBayesLM
library(R.filesets)
## Utils 
## Summarized experiment 

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

safeBPParam <- function(nworkers=14) {
    if (.Platform$OS.type=="windows") {
        BiocParallel::SerialParam()
         BiocParallel::SnowParam(nworkers)
    } else {
        BiocParallel::MulticoreParam(nworkers)
    }
}

### map event to months 

EVENT_MAP=list('SC' = -3,  'BL' =  0,  'V01'=3,    'V02'=6,    'V03'=9,    'V04'=12,   'V05'=18,   'V06'=24,   'V07'=30,   
               'V08'=36,    'V09'=42,    'V10'=48,    'V11'=54,   'V12'=60,   'V13'=72,   'V14'=84,   'V15'=96, 'V16'=108, 'V17'=120,'V18'=132,'V19'=144 , 
               'ST' = -6)


#diff_vars<-colnames(samples_metadata(MOFAobject))[grep('diff', colnames(samples_metadata(MOFAobject)))]
#diff_vars<-diff_vars[!grepl('clust',diff_vars)]
#diff_vars
selected_covars_broad<-c('COHORT', 'AGE', 'SEX','NP1RTOT', 'NP2PTOT','NP3TOT', 'updrs3_score_on', 
                         'NP1_TOT', 'NP2_TOT','NP3_TOT', 'NP4_TOT',
                         'NHY', 'NP3BRADY',
                         'NP3RIGN', 'SCAU5', 'MCATOT','moca', 
                         'MCAVFNUM', 'MCACLCKH', 'cogstate','sft' , 'VLTFRUIT', 
                         'ptau', 'asyn', 'tau_ab', 'tau_asyn', 'abeta', 'ess_cat', 
                         'HVLTRDLY',
                         'PDSTATE', 'NP3RTCON', 
                         'stai_state', 'stai_trait'  ,'STAIAD26', 'NP1ANXS', 'NP3GAIT', 
                         'SCAU7', 'NP3SPCH', 'NP3RISNG', 'NP2EAT', 
                         'NP3RTARU', 'RBD_TOT', 
                         'con_putamen', 
                         'td_pigd_old_on', 'PD_MED_USE' , 'Outcome', 
                         'rigidity','months', 
                         'con_putamen', 'con_putamen_V10', 
                         'change', 
                         'asyn', 'CSFSAA', 'NP3_TOT_LOG_SCALED', 
                         'NP3_TOT_diff_V16', 'SCAU_TOT_diff_V16', 'NP2_TOT_diff_V16',
                         'con_putamen_diff_V10', 'hi_putamen_diff_V10',
                         'MCA_TOT_diff_V16', 'SITE', 'Plate','Usable_bases_SCALE', 
                         'Neutrophils....', 'Lymphocytes....', 'Neutrophils.Lymphocytes', 'RIN', 
                         'Uniquely.mapped')
                   #      diff_vars)
#'DYSKIRAT')


# sPLIT DIAGNOSIS vs progression  
selected_covars2<-c( 'AGE', 'SEX',
                     'NP2PTOT','NP3TOT',
                     'updrs2_score','updrs3_score', # Todo some are missing from these scores...
                    'NHY', 
                     'NP3RIGN',
                     'rigidity', 
                     'NP3RTARU',
                     'PDSTATE',  
                     'td_pigd_old_on', 
                     'PD_MED_USE' , 
                     'months', 
                     'con_putamen_V10', 
                     'CSFSAA', 
                     'mean_striatum_V10', 
                     'ab_asyn')

selected_covars2_progression<-c( 'AGE', 'SEX',
                                 #'NP1_TOT', 
                                 'NP2PTOT_LOG','NP3TOT_LOG',
                                 'updrs3_score',
                                 #'NP4_TOT',
                                 'NHY', 
                                 # 'SCAU5',
                                 # NON-MOTOR
                                 'moca','sft', 'VLTFRUIT', 'VLTVEG', 'PUTAMEN_R_V10', 
                                 # CSF BIOMARKERS 
                                 'abeta_LLOD', 'HVLTRDLY',
                                 #  'scopa', 
                                 # 'stai_state', 'stai_trait', 
                                 'rigidity', 
                                 'NP3RTARU',
                                 # not significant: 
                                 #  'ptau',    'ab_asyn', 
                                 'PDSTATE',  
                                                                  'PD_MED_USE' , 
                             #    'months', 
                                 'con_putamen_V10', 
                             'asyn' , 'CSFSAA', 
                                 'mean_striatum_V10', 
                                 
                                 ## WHICH factors have to do with the change in the scale
                                 # And the change in the datascan binding  in the future?
                                 # THESE factors are the ones that we actually WANT                                #  'con_putamen_diff_V10', 'hi_putamen_diff_V10',
                                 'SITE', 'Plate', 
                                 'Neutrophil.Score', 
                                 
                                'RIN.Value', 'Multimapped....', 'Uniquely.mapped....'
                                 
                                 #'MCA_TOT_diff_V16', 
                                 
                                 #'NP3_TOT_LOG_SCALED', 
                                 #'RBD_TOT_diff_V16', 'MCATOT_diff_V16'
)



mt_kv<-read.csv(paste0(output_files, 'metadata_key_value.csv'), header = FALSE)
mt_kv$V1<-gsub(' |\'|\"','',mt_kv$V1 )
mt_kv$V2<-gsub(' |\'|\"','',mt_kv$V2 )


## GENES RELATED TO THE BATCH EFFECT
filtered_genes<-read.csv(paste0(data_dir, 'ppmi/ppmi_data/rnaseq/filteredGenes.csv'))

remaining<-read.csv('ppmi/remaining_genes.csv', header = FALSE)
remaining_genes<-gsub( '\\..*', '', remaining$V1)

batch_effect_genes<-filtered_genes$perc0.1
remove_genes<-batch_effect_genes

load_se_all_visits<-function(input_file, combined){
  #'
  #' @param input_file file cntaining rnas or mirnas counts for all visits 
  #' 
  #' 
  #write.csv2(rna_all_visits, gz1, row.names = TRUE)
  
  raw_counts<-as.matrix(fread(input_file, header=TRUE), rownames=1)

  
  raw_counts_all<-raw_counts
  class(raw_counts_all) <- "numeric"
  ## They seem to have taken averages for replicas so need to fix
  raw_counts_all<-round(raw_counts_all)
  ## Question: why are there duplicate samples - seems to be controls!
  ## first filter what is in metadata and mirnas ?
  
 # grep('1009',colnames(raw_counts_all))
  
  
  se<-getSummarizedExperimentFromAllVisits(raw_counts_all, combined)
  return(se)
}

### TODO: move to a utils / preprocessing file because it is used also for proteoomics
library(SummarizedExperiment)


getSummarizedExperimentFromAllVisits<-function(raw_counts_all, combined){
  #'
  #'
  #' key id includes patno_event_id + ON_OFF 
  #'
  #'
  
  
  raw_counts_all<-raw_counts_all[,!duplicated(colnames(raw_counts_all), fromLast=TRUE)]
  combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)
  ### some samples do not exist in metadata so filter them out 
  ## 
  common_samples<-intersect(colnames(raw_counts_all),combined$PATNO_EVENT_ID)
  unique_s<-colnames(raw_counts_all)[!(colnames(raw_counts_all) %in% common_samples)]
  # TODO: replace with function: fetch_metadata_by_patient_visit- test first
  metadata_filt<-fetch_metadata_by_patient_visit(common_samples)
  #metadata_filt<-combined[match(common_samples, combined$PATNO_EVENT_ID),]
  raw_counts_filt<-raw_counts_all[,match(common_samples, colnames(raw_counts_all))]
  dim(metadata_filt)[1] ==dim(raw_counts_filt)[2]

  
  #subset sample names

  se=SummarizedExperiment(raw_counts_filt, colData = metadata_filt)
  
  metadata_filt$COHORT_DEFINITION
  
  return(se)
  
  
}


#se_filt<-se_filt_all[[cluster_id]]
#dim(assay(se_filt))
#se_filt$COHORT
preprocess_se_deseq2<-function(se_filt, min.count=10){
  #' 
  #' Preprocess metadata of summarized experiment 
  #' 
  #'  
  
  # Raw data preprocess 
  # 
 
  ## Turn to factors for deseq
  
  se_filt$SITE <- as.factor(se_filt$SITE); dim(assay(se_filt))
  se_filt$Neutrophils<-scale(se_filt$`Neutrophils....`)
  se_filt$Lymphocytes<-scale(se_filt$`Lymphocytes....`)
  se_filt$Neutrophil.Score<-scale(se_filt$Neutrophil.Score)
  
  
  se_filt$Usable_Bases_SCALE<-as.numeric(scale(se_filt$Usable.Bases....))
  se_filt$SEX<-as.factor(se_filt$SEX)
  se_filt$SITE<-as.factor(se_filt$SITE)
  se_filt$Plate<-as.factor(se_filt$Plate); dim(assay(se_filt))
  
  if ('Neutrophils.LD' %in% colnames(colData(se_filt))){
    colData(se_filt)[, colnames(estimations)]<-scale(colData(se_filt)[, colnames(estimations)])
  }



  #se_filt$COHORT[ which(is.na(se_filt$COHORT))]<-'Unknown'
  se_filt<-se_filt[ , !is.na(se_filt$Neutrophil.Score)]

  se_filt<-se_filt[ , !is.na(se_filt$COHORT)]
  
  se_filt<-se_filt[ , !is.na(se_filt$Usable_Bases_SCALE)]
   dim(assay(se_filt))
  
  ### OUTPUT THE FILTERED se_filt 
  ind<-which(is.na(se_filt$AGE_AT_VISIT))
  se_filt[,ind]$AGE_AT_VISIT<-get_age_at_visit(colData(se_filt[,ind]))
  se_filt$AGE_AT_VISIT<-scale(se_filt$AGE_AT_VISIT)
  se_filt$AGE_SCALED[is.na(se_filt$AGE_SCALED)]<-mean(se_filt$AGE_SCALED, na.rm=TRUE)
  se_filt<-se_filt[,!(is.na(se_filt$SEX))]
  
  table(colData(se_filt)[,c( 'EVENT_ID', 'SEX')])
  
  colData(se_filt)[,c( 'EVENT_ID', 'SEX', 'AGE', 'PATNO')]
  
  
  #### remove genes
  #se_filt<-se_filt[!(rownames(se_filt) %in% remove_genes),]
  # removed genes associated with the batch effects as explained in Craig 2020 paper
  rownames(assay(se_filt))
  
  if (!process_mirnas){
    # filter by specified genes in
      assay_r<-gsub( '\\..*','' ,rownames(assay(se_filt)))
      se_filt<-se_filt[assay_r %in% intersect(assay_r, remaining_genes),]

  }
 
  
  # 
  se_filt<-filter_se_byExpr(se_filt, min.count=min.count);   dim(assay(se_filt) )


 
  min.count
 
  
  return(se_filt)
}



deseq_by_group<-function(se_filt, formula_deseq, min.count=10,cell_corr_deseq=TRUE){
  #'
  #' @param 
  # TODO: add plate and/OR site 
  # se_filt1 neutrophil counts, and usable bases
  
  #IF NEUTROPHILS IN DESIGN
  if (cell_corr_deseq){
    #se_filt<-se_filt[,!(is.na(se_filt$`Neutrophils....`))] ;dim(assay(se_filt))
    se_filt<-se_filt[,!(is.na(se_filt$Neutrophil.Score))] ;dim(assay(se_filt))
    
    #se_filt<-se_filt[,!(is.na(se_filt$`Lymphocytes....`))] 
  }
  # If using kmeans grouping filter if NA
  
  
  se_filt<-se_filt[,!(is.na(se_filt$COHORT))] 
  
  # preprocess data turn to factors etc
  se_filt<-preprocess_se_deseq2(se_filt, min.count=min.count)
  assay(se_filt)
  # if correcting by medication 
  se_filt$PDMEDYN = as.factor(se_filt$PDMEDYN)
  se_filt$PDMEDYN[is.na(se_filt$PDMEDYN)]=0
  se_filt$LEDD[is.na(se_filt$LEDD)]=0 # add zeros to na then scale!
  se_filt$LEDD_scaled<- scale(se_filt$LEDD)
  
  
  assay(se_filt)
  
  
  ddsSE <- DESeqDataSet(se_filt, 
                        design = as.formula(formula_deseq))
  ddsSE<-estimateSizeFactors(ddsSE)
  
  #vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
  
  deseq2Data <- DESeq(ddsSE, parallel=TRUE, BPPARAM = safeBPParam())
  #deseq2Data <- DESeq(ddsSE, parallel=TRUE)
  #deseq2Data <- DESeq(ddsSE)
  
  deseq2Results<-results(deseq2Data)
  deseq2ResDF <- as.data.frame(deseq2Results)
  
  padj_T_hv<-0.05
  log2fol_T_hv<-0.1
  
  ### this is also done later on -- save from there? 
  deseq2ResDF$mofa_sign<- ifelse(deseq2ResDF$padj <padj_T_hv & abs(deseq2ResDF$log2FoldChange) >log2fol_T_hv , "Significant", NA)
  deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
  return(deseq2ResDF)
}


#remove_genes


fetch_metadata_by_patient_visit<-function(patno_event_ids, combined=combined_bl_log, PDSTATE=NULL){
  

   #'
   #' @param PATNO_EVENT_ID
   #'
   #' fetch one row per patient back -- if multiple NP3TOT it brings the maximum! 
   #' 
   #' 
   # TODO: add more criteria to filter eg. pdstate='OFF'; pag_name_m3='NUPDRS3'
   # PRIORITIZE OFF ?  - KEEP ALSO ON FOR THE MISSING ONES?? SO WE CAN HAVE THE OTHER DATA? 
   
   #combined_bl_log_sel<-combined_bl_log_sel %>%
  #   group_by(NP3_TOT) %>%
  #   summarize(across(everything(), max))%>%
  #   as.data.frame()
  
  # patno_event_ids = paste0(patnos, event_ids)
   #combined<-combined_bl_log

   
   
    metadata_filt_dups<-combined[combined$PATNO_EVENT_ID %in% patno_event_ids ,]
 
    metadata_filt_dups<-as.data.table(metadata_filt_dups)
    
    max_np3_unique<-metadata_filt_dups[metadata_filt_dups[, .I[which.max(NP3TOT)], by=PATNO_EVENT_ID]$V1]
    
    
    dim(max_np3_unique)
    length(patno_event_ids)

    #missing<-metadata_filt_dups[is.na(metadata_filt_dups$EVENT_ID)]
    #max_np3_unique<-rbind(max_np3_unique, missing)
    # TODO: mca tot is zero.. 
      # as.data.frame(metadata_filt_dups[, c('PATNO_EVENT_ID','MCA')])
    

    dups<-duplicated(metadata_filt_dups[, c('PATNO_EVENT_ID', 'PDSTATE', 'NUPSOURC', 'PAG_NAME_M3', 'PAG_NAME_M4','NP3TOT')]%>%
     arrange(PATNO_EVENT_ID, PDSTATE, NUPSOURC, PAG_NAME_M3, NP3TOT))
    
    metadata_filt_dups<-metadata_filt_dups[!dups,]
    
    metadata_filt_dups[, c('PATNO_EVENT_ID', 'PDSTATE', 'NUPSOURC', 'PAG_NAME_M3', 'PAG_NAME_M4','NP3TOT')]%>%
                       arrange(PATNO_EVENT_ID)
    
    dups_strict<-duplicated(metadata_filt_dups[, c('PATNO_EVENT_ID')])
    metadata_filt_unique<-metadata_filt_dups[!dups_strict,]
     
    
    ### if we cannot find it for all then use what was given? 
    # match to ensure same order
    max_np3_unique<-as.data.frame(max_np3_unique[match(patno_event_ids, max_np3_unique$PATNO_EVENT_ID ), ])
    dim(max_np3_unique)
    
    #max_np3_unique[,'PATNO_EVENT_ID']<-patno_event_ids
    ### HERE IF IT DID NOT FIND THE metadata it adds back the id 
    ### replace the ones with missing NP3TOT
    # TODO: fix here what to do for the ones without np3 total--> sum to obtain it? 
    missing<-patno_event_ids[is.na(max_np3_unique$EVENT_ID) ]
    missing_metadata<-as.data.frame(combined[match(missing, combined$PATNO_EVENT_ID ), ])
    max_np3_unique[is.na(max_np3_unique$EVENT_ID),]<-missing_metadata

    missing2<-patno_event_ids[is.na(max_np3_unique$EVENT_ID) ]

    
  return(max_np3_unique)

 }
 

 ##########
 ######

 imaging_variables_diff<-c('updrs3_score', 'updrs3_score_on', 'updrs3_score_LOG', 'updrs3_score_on_LOG', 'con_putamen', 'hi_putamen',
                           'updrs2_score', 'updrs2_score_LOG', 'moca' , 'scopa')
 scale_vars_diff=c('NP3TOT', 'NP2PTOT', 'RBD_TOT', 'MCATOT' ,'SCAU_TOT' )### todo add upsit and other scales? 
 
 
 get_diff_zscores<-function(patno_event_ids,imaging_variables_diff,scale_vars_diff ){
            #### obtains the zscore of the changes for the specific group supplied
              ### MEAN AND sd for scaling depends on the group!
   #' 
   #' 
   #' @param patno_event_ids
   #' could also supply the specific variables to diff

   
   # 

            df_all<-fetch_metadata_by_patient_visit(patno_event_ids , combined=combined_bl_log)

             t1<-'BL';  t2='V10';
             df_all$moca_V10
             #df=df_all
             paste0(imaging_variables_diff, '_BL')%in% colnames(df_all)
             df_change1= get_changes(df_all,imaging_variables_diff, t1, t2 )
             #df_all<-cbind(df_all, df_change1)
            
             t1<-'BL';  t2='V12';
             df_all$moca_V10
             #df=df_all
             df_change_V12= get_changes(df_all,imaging_variables_diff, t1, t2 )
             
             
             
              t1<-'BL';  t2='V16';
             colnames(df_all)[grep('V16', colnames(df_all))]
             ### hack here warning!! I added the data from the curated because it was missing 
             df_all$MCATOT_BL <- df_all$moca_BL
             
             
             # TODO: for mca tot it is not available so add moca for baseline and MCATOT for other visits 
             
             not_in_df<-scale_vars_diff[!scale_vars_diff %in% colnames(df_all)]
             if (length(not_in_df)>0){
               print(paste('ERROR: missing scales: ', not_in_df))
             }
             df_change2= get_changes(df_all,scale_vars_diff, t1, t2 )
             
             df_all[,paste0(scale_vars_diff, '_V14')]
             paste0(scale_vars_diff, '_V14') %in% colnames(df_all)
             df_change_V14= get_changes(df_all,scale_vars_diff, t1='BL', t2='V14' )
             
             scale_vars_diff
             
             df_change_av=get_av_change(df_all,scale_vars_diff,'BL','V13', 'V14' )
             
             
             
            df_change_total<-cbind(df_change1, df_change2,df_change_V12,df_change_V14, df_change_av)
              
        return(df_change_total)
        
        
 }

## Create the summarized experiment by selecting VISITS and cohorts 
filter_se<-function(se, VISIT, sel_coh, sel_sub_coh=FALSE){
  
  #' Takes the raw file with all counts
  #' Filters summarized experiment by selecting VISITS and cohorts 
  #' @param VISIT
  #' @param sel_coh
  
  ##### 2.   start filtering the experiment  to normalize as appropriate 
  ## Option 1: normalize cohort and EVENT separately!! 
  # ALSO MAKE SURE THAT they are in cohort in the conversion cohort too!!
  
  
  if (length(sel_subcoh)==1 && sel_subcoh==FALSE){
        se_filt<-se[,((se$EVENT_ID %in% VISIT) & (se$COHORT %in% sel_coh ) & (se$CONCOHORT %in% sel_coh ))]
    
  }else{
      se_visit<-se[,se$EVENT_ID %in% VISIT]

      if (1 %in% sel_coh){ 
        ids<-c(se_visit$INEXPAGE %in% sel_subcoh ) ## filter the ids in parkinsons 
      }
      if (2 %in% sel_coh){
        ids<- c(ids | se_visit$INEXPAGE %in% 'INEXHC') ## also extract controls 
        
      }
      se_filt<-se_visit[, ids]
      
      
    
  }
     
  dim(se_filt)     

  
  Sample<-colnames(se_filt)
  sample_info<-DataFrame(Sample=Sample)
  
 
  
  ##### Define
  
  ### TODO: Question: Should I input everything into the matrix to normalize? 
  ### And then filter 
  
  ### batch effect and normalization 
  # Create a separate matrix with counts only
  # Include batch information if there is any
  #sample_info$Batch <- as.factor(sample_info$Batch)
  
  
  ### DEFINE THE DESEQ OBJECT with the groups appropriately 
  se_filt$EVENT_ID=as.factor(se_filt$EVENT_ID)
  se_filt$COHORT=as.factor(se_filt$COHORT)
  se_filt$PATNO=as.factor(se_filt$PATNO)
  
  
  return(se_filt)
  
}




### METADTA 

#install.packages('eeptools' )
#library('eeptools')

#as.Date(as.character(new$STATUS_DATE), format = "MM/YY",)

# TODO: fix 
get_age_at_visit<-function(new){
  #'
  #'
  #'
  AGE_AT_VISIT<-as.numeric(gsub('.*/','',new$STATUS_DATE)) - as.numeric(gsub('.*/','',new$BIRTHDT))
  return(AGE_AT_VISIT)
  }
#x_age <- age_calc( as.Date(new$BIRTHDT),          # Convert birth to age
##                   as.Date(new$STATUS_DATE),
#                  units = "years")


#### data specifc 
get_symbols_vector<-function(ens ){
  #' @param ens ensemble ids to conver to symbols 
  #' @returns symbols_ordered the total 
  #'  
  #'  
  ens<-gsub('\\..*','', ens ) # remove if there is something after the dot 

  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  symbols <- symbols[!is.na(symbols)]
  symbols_ordered <- symbols[match(ens, names(symbols))]
  na_ind<-is.na(symbols_ordered);
  
  # Add ensembl ids if no symbol found
  symbols_ordered[na_ind]=ens[na_ind]
  return(symbols_ordered)
  
  
}



######## DE ANALYSIS #######
#results_de<-mark_signficant(
  # test
 # de_res= deseq2ResDF
##deseq2ResDF$padj
  #padj_T = padj_T_overall; log2fol_T = log2fol_T_overall; padj_name ='padj'
  #log2fc_name = 'log2FoldChange'   
  #outdir_single=outdir_s_p

mark_significant<-function(de_res, padj_T, log2fol_T, padj_name='padj', log2fc_name='log2FoldChange', outdir_single=outdir_s ){
  ## mark a significant column and write to file
  
  signif_file<-paste0(outdir_single,'/significant', padj_T, '_',log2fol_T, '.csv')
  
  de_res$significant <- ifelse(de_res[,padj_name] < padj_T , "Significant", NA)
  de_res$sign_lfc <- ifelse(de_res[,padj_name] < padj_T & abs(de_res[,log2fc_name]) >log2fol_T , "Significant", NA)
  # Examine this data frame
  # Order the significant to save as a new output file 
  head(de_res)
  which(de_res$significant=='Significant')
  # LARGER ONE not saved 
  sign_only<-de_res[de_res$sign_lfc=='Significant',]
  sign_only_ordered<-sign_only[order(sign_only[,padj_name], decreasing = FALSE),]
  sign_only_ordered<-sign_only[order(-sign_only[,'abslog2pval'], decreasing = FALSE),]
  
  sign_only_ordered<-sign_only_ordered[!is.na(sign_only_ordered$sign_lfc=='Significant'),]

  which(sign_only_ordered$significant=='Significant')
  write.csv(sign_only_ordered,signif_file, row.names = TRUE)
  ### create also a more strict file? 
  return(de_res)
}



######## ENRICHMENT ANALYSIS 


get_ordered_gene_list<-function(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 ){
  
  #### Gives a gene list cut by the thresholds padj_T and log2fol and orders by specific metric supplied 
  #' @param padj_T filter the genes by metric padj_T
  #' @param log2fol_T  filter the genes by metric  log2fol_T, default: 0
  #' @param order_by_metric metric to order the gene list by 
  res=deseq2ResDF
  #res=deseq2ResDF_2
  
  res$sign_lfc <- ifelse(res$padj <padj_T & abs(res$log2FoldChange) >log2fol_T , "Significant", NA)
  
  length(which(!is.na(res$sign_lfc )))
  res=res[res$sign_lfc=='Significant'& !is.na(res$sign_lfc),]
  
  
  # Order the DE gene list by the stat statistic 
  #remove negatives thatw ere introduced with vst transofrmations
  
  res$log2pval<-res$log2FoldChange*-log10(res$padj)
  res$signlog2pval<-sign(res$log2FoldChange)*-log10(res$padj)
  res<-res[res$baseMean>0,]
  
  #res <- res[order(-res$stat),]
  res <- res[order(-res[,order_by_metric]),]
  gene_list<-res[, order_by_metric]
  names(gene_list)<-rownames(res)
  
  return(gene_list)
  
  
  
}


#### enrichment packages 
#library('GOfuncR')
require(DOSE)
library(clusterProfiler)
library(ggplot2)

suppressWarnings(library('R.filesets' ))
library('enrichplot' )
require(DOSE)
library('enrichplot')


#deseq2ResDF
#results_file_cluster
#N_DOT=20, N_EMAP=30 , N_NET=10)

run_enrich_per_cluster<-function(deseq2ResDF, results_file,N_DOT=15, N_EMAP=25, N_NET=15){
  #'
  #' Get the ordered list 
  #' @param deseq2ResDF
  #' 
  #'
  #deseq2ResDF=deseq2ResDF1
  gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )
  gene_list_ord=gene_list
  names(gene_list)<-gsub('\\..*', '',names(gene_list))
  gse=run_enrich_gene_list(gene_list, results_file)
 
  return(gse)
}





run_enrich_gene_list<-function(gene_list, results_file, N_DOT=15,N_EMAP=30, pvalueCutoff_sig=0.05, pvalueCutoff=1){
  #'
  #' Run enrichment GSEA using an ordered gene list, write, and plot and save the results 
  #' @param gene_list
  #' ordered gene list to run gse 
  #' @return return the significant results 
  #' 
  #'

  #pvalueCutoff=1; gse_full; N_DOT=15;results_file=results_file_cluster;
    #results_file


  gse_full <- clusterProfiler::gseGO(gene_list, 
                                     ont=ONT, 
                                     keyType = 'ENSEMBL', 
                                     OrgDb = 'org.Hs.eg.db', 
                                     pvalueCutoff  = pvalueCutoff  )
  
  
  gse = write_filter_gse_results(gse_full, results_file, pvalueCutoff)
  
  gse = dplyr::filter(gse_full, p.adjust < pvalueCutoff_sig)
    if  (dim(gse)[1]>2 ){
      # check if the enrichment returned more than one resulting paths before plotting!

      enrich_plots <- run_enrichment_plots(gse=gse,results_file=results_file, N_DOT=N_DOT, N_EMAP=N_EMAP, N_NET=15,
          run_ORA=FALSE, geneList = gene_list)
    }
  
  return(gse)
}



#results_file_ora = paste0(results_file_cluster,'_ora' )

run_ora_gene_list<-function(gene_list_ord, results_file_ora, keyType='SYMBOL', N_DOT=15, N_EMAP=25){
  #'
  #' Run enrichment and write the results 
  #' @param gene_list
  #' ordered gene list to run gse 
  #'
  #' 
  #'
  #gene_list_ord=gene_list
  #keyType='ENSEMBL'
  gene_list_ord_abs=abs(gene_list_ord)
  
  
  gene_list_ord_abs_ora=gene_list_ord_abs[order(gene_list_ord_abs, decreasing = TRUE)]
  hist(gene_list_ord_abs_ora)
  high_quant<-quantile(gene_list_ord_abs_ora, 0)
  high_quant_no<-length(gene_list_ord_abs_ora[gene_list_ord_abs_ora>high_quant])
  
  gse_protein_full_enrich <- clusterProfiler::enrichGO(names(gene_list_ord_abs_ora[1:80]), 
                                                       ont=ONT, 
                                                       keyType = keyType, 
                                                       OrgDb = 'org.Hs.eg.db', 
                                                       pvalueCutoff  = pvalueCutoff,
                                                       pool =FALSE)
  
  
  
  
  ### to run mofa results
  
  process_mirnas=FALSE
  process_mofa=TRUE
  run_ORA=TRUE
  
  gse_full=gse_protein_full_enrich
  pvalueCutoff_sig = 0.1
  gse=write_filter_gse_results(gse_full, results_file_ora, pvalueCutoff=0.1)
  write.csv(as.data.frame(gse_full@result), paste0(results_file_ora, '.csv'))
  
  gse=dplyr::filter(gse_full, p.adjust < pvalueCutoff_sig)
  gse@result$Description
  
  gse@result$p.adjust
  
  #N_DOT=15
  #N_EMAP=30
  if  (dim(gse)[1]>2 ){
    
    
       print(paste(factor,'sig'))
         which(gse_full@result$p.adjust<0.05)
  
        enrich_plots<-run_enrichment_plots(gse=gse,results_file=results_file_ora, N_DOT=N_DOT, N_EMAP=N_EMAP, run_ORA=TRUE)
  }
  return(gse_full)
}




write_filter_gse_results<-function(gse_full,results_file,pvalueCutoff  ){
  
  ### Takes the full gse results, ie. without threshold significance, 
  # saves it, 
  # filters it by pvalueCutoff_sig
  # and saves the filter 
  #' @param gse_full full gse results objects 
  #' @param results_file the file name to write results  (without .csv)
  #' @param pvalueCutoff the pvalue used to obtain the gse results 
  
  write.csv(as.data.frame(gse_full@result), paste0(results_file, pvalueCutoff, '.csv'))
  pvalueCutoff_sig<-0.05
  gse_sig_result<-gse_full@result[gse_full@result$pvalue<pvalueCutoff_sig,]
  write.csv(as.data.frame(gse_sig_result), paste0(results_file, pvalueCutoff_sig, '.csv'))
  
  # rewrite
  dim(gse_full); dim(gse_sig_result)
  ## filter gse result to significant only 
  gse=dplyr::filter(gse_full, p.adjust < pvalueCutoff_sig)
  return(gse)
}


plot_enrich_compare<-function(gse_compare,enrich_compare_path, N_EMAP=35, N_DOT=8){
  ### GSE COMPARE ANALYSIS 
  dot_comp<-clusterProfiler::dotplot(gse_compare, showCategory=N_DOT, split=".sign") + facet_grid(.~.sign)
  dot_comp
  ggsave(paste0(enrich_compare_path, 'dt','.jpeg' ), plot=dot_comp,
         dpi=250, width=12, height=10, 
  )
  
  # N_EMAP 
  #N_EMAP=80
  gse_compare_x <- enrichplot::pairwise_termsim(gse_compare)
  
  emap_comp<-emapplot(gse_compare_x, showCategory=N_EMAP,
                      cex.params = list(category_label = 0.9) ) 
  emap_comp
  ggsave(paste0(enrich_compare_path, 'emap',N_EMAP,'.jpeg' ), plot=emap_comp,
         dpi=200, width=10, height=10, units='in')
  


## first convert gene ID to Symbol
gse_compare.orig <- setReadable(gse_compare, 'org.Hs.eg.db', 'ENSEMBL')

## use new way of specifying visualization options
cex.params = list(category_label = 0.6, gene_label = 0.4)
  
cnet_comp<-cnetplot(gse_compare.orig, showCategory = 10, 
         cex.params = cex.params)

ggsave(paste0(enrich_compare_path, 'cnet.jpeg' ), plot=cnet_comp,
         dpi=200)
  
  
}

#gse=gse_mofa_rna
#results_file = results_file_mofa
#N_EMAP=50
#geneList =NULL  
#title_p

#gse=gse_mofa_sig
#gse=gse1
#N_EMAP=25; N_DOT=15; N_TREE=16; N_NET=30
#results_file = results_file_mofa, N_EMAP=50,geneList =NULL  )
library('clusterProfiler')
#options(warn=0, error=NULL)
run_enrichment_plots<-function(gse, results_file,N_EMAP=25, N_DOT=15, N_TREE=16, N_NET=20, showCategory_list=FALSE,
                               process_mofa=FALSE, text_p='', title_p='', geneList=NULL, run_ORA=FALSE){

  #gse=gse_mofa_rna; 
  #geneList=gene_list_ord_g
  #N_EMAP=25; N_DOT=15; N_TREE=16; N_NET=30
  #run_ORA=FALSE
  N=25
 
    write_n=FALSE
  
  
  ### print a signed and unsigned dotplot 
  # because it does not make sense if we dont rank by logFC
  # or in the mofa case where we rank by importance in factor 
   # N_DOT=100
 
  dp<-clusterProfiler::dotplot(gse, showCategory=N_DOT, 
              font.size=15
  )
  dp
  dp<-dp+theme(axis.ticks=element_blank() , 
               axis.text.x = element_blank(),
               plot.caption= element_text(hjust = 0, face = "italic", size=20), 
               plot.title=element_text(size=20)) +
    labs(caption=text_p, title=title_p)

  
    
  
  if (isRStudio){
    show(dp)
    
  }
  
  
  
  if (process_mirnas){
    width=6}else{width=6}
  
  ggsave(paste0(results_file, '_dot', N_DOT, '.jpeg'), 
         plot=dp, width=width, height=N_DOT*0.5, 
         dpi = 200)
  
  if (!(process_mirnas) & !(run_ORA)){
    
      print('signed dotplot')
      #dp_sign<-  clusterProfiler::dotplot(gse, showCategory=N_DOT, split=".sign") + facet_grid(.~.sign)
      #ggsave(paste0(results_file, '_dsign', N_DOT,  '.jpeg'), width=8, height=N*0.7)
    
   
    
  }

  #### EMAP PLOT 
  options(ggrepel.max.overlaps = Inf)
  x2 <- pairwise_termsim(gse )
  #if (process_mirnas){N=15}
  p<-emapplot(x2,showCategory = N_EMAP,
              layout.params = list(layout = "nicely"), 
              cex_label_category=0.8, 
                )
  p_enrich <- p + theme(text=element_text(size=12))
  p_enrich
  
  if (is.numeric(N_EMAP)){write_n=N_EMAP}
  ggsave(paste0(results_file, '_emap_', write_n,  '.jpeg'), width=9, height=9, 
         dpi = 200)
  
  
  #### Ridge plot: NES shows what is at the bottom of the list
  
  
  N_RIDGE=25
  # only if all 3 are false run it 
  if ( !(process_mirnas) & !(process_mofa) & !(run_ORA)){
    print('ridge')
    r_p<-ridgeplot(gse, showCategory = N_RIDGE)
    if (dim(r_p$data)[1]>0){
      r_p
      ggsave(paste0(results_file, '_ridge_', N_RIDGE, '.jpeg'), width=8, height=8)
      
    }
  }
  
  
  
  
  #### Gene-concept plot 
  
  
  
  if (gse@keytype=='ENSEMBL' ){
    gse_x <- setReadable(gse, 'org.Hs.eg.db', 'ENSEMBL')
    
  }else{
    gse_x=gse
    
  }
  if (!is.null(geneList)){
#  geneList=gene_list1
  
    node_label<-"gene"
    node_label<-"category"
    node_label<-"all"
    graphics.off()
    
    p2_net<- cnetplot(gse_x,
                      node_label=node_label,
                      cex_label_category = 0.9,
                      showCategory=N_NET, 
                      foldChange =  geneList)
    
    p2_net
    show(p2_net)
    #graphics.off()
    if (is.numeric(N_NET)){write_n=N_NET}
    unlist(strsplit('all', '',1))[[1]]

    ggsave(paste0(results_file, '_gc_', unlist(strsplit(node_label, ''))[[1]], '_',write_n, '.jpeg'), width=20, height=20)
  }else{
    p1_net <- cnetplot(gse_x)
    show(p1_net)
    ggsave(paste0(results_file, 'gc', '.jpeg'), width=12, height=12)

    
  }
  
  
  node_label='all'
  as.character(node_label)
  
  ####Visualize go terms as an undirected acyclic graph 0
  #if (!process_mirnas){
  #  goplot(gse_x)
  #  ggsave(paste0(results_file, '_goplot_', node_label, '_',write_n, '.jpeg'), width=8, height=8)
    
  #} 
  library(ggtree)
  library(ggplot2)
  
  #install.packages('ggtree')
  
  #### heatmap
  plot_tree=FALSE
  if (plot_tree){
    nCluster=ifelse(dim(x2)[1]<4,1, 4) 
    p1 <- treeplot(x2,showCategory =N_TREE, nWords=0, nCluster=nCluster)
    
    
    
    p2_tree <- treeplot(x2, hclust_method = "average", 
                        showCategory =N_TREE, nWords=0, nCluster=nCluster,
                        label_format =50, 
                        fontsize = 300, 
                        extend=-0.001, 
                        offset=15, 
                        hilight=FALSE, 
                        branch.length=0.1)
    
    #aplot::plot_list(p1, p2_tree, tag_levels='A')
    #ggsave(paste0(results_file, '_clusterplot_', node_label, '_',N, '.jpeg'), width=8, height=8)
    
    # p2_tree<-p2_tree+theme(plot.caption= element_text(hjust = 0, face = "italic", size=20)) +
    # labs(caption=text_p, title=title_p)
    #write_n='test'
    #N_TREE=10
    ggsave(paste0(results_file, '_cp_',N_TREE, '.jpeg'),
           width=10, height=0.5*log(N_TREE)+3, dpi=200)
    #
    
  }
  
      
    
  
  return(list(dp, p_enrich))
  
}






get_genelist_byVisit<-function(VISIT){
  
  ### Input visit AND return list 
  ##'
  ##'
  
  deseq2ResDF_2 = as.data.frame(read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1))
  gene_list<-get_ordered_gene_list(deseq2ResDF_2,  order_by_metric, padj_T=1, log2fol_T=0 )
  names(gene_list)<-gsub('\\..*', '',names(gene_list))
  return(gene_list)
}



#### mirna enrich ####


#install.packages('VennDiagram')
library(rbioapi)
library('VennDiagram')
library(enrichplot)




#table(mieaa_all_gsea$Category)
Category<-'Gene Ontology (miRWalk)';
Category<-'Annotation (Gene Ontology)';
Category<-'GO Biological process (miRPathDB)';




#install.packages("remotes")
#remotes::install_github("jmw86069/jamenrich")
library('multienrichjam')
library('clusterProfiler')

mirna_enrich_res_postprocessing=function(mieaa_all_gsea,mir_results_file,  Category='GO Biological process (miRPathDB)'){
#mirna_enrich_res_postprocessing=function(mieaa_all_gsea, Category='GO Biological process (miRPathDB)',mir_results_file){

  #' post-process mieaa enrichment analysis results 
  #' convert to enrich result to be able to use with cluster profiler plotting functions
  #' @param mieaa_all_gsea output from mieaa
  #' @param Category which category to choose from enrich dbs eg.  'GO Biological process (miRPathDB)'
  #'
  colnames(mieaa_all_gsea)<-gsub('-','.', colnames(mieaa_all_gsea))
  colnames(mieaa_all_gsea)<-gsub('/','.', colnames(mieaa_all_gsea))
  
  
  
  mieaa_gsea_1<-mieaa_all_gsea[mieaa_all_gsea$Category==Category,]
  ## write output results 
  write.csv(mieaa_gsea_1, paste0(mir_results_file, '_', pvalueCutoff,'.csv' ))
  
  ### cut filtered datasets only
  mieaa_gsea_1_cut<-mieaa_gsea_1[mieaa_gsea_1$P.adjusted<Padj_T_paths, ]
  mir_paths<-mieaa_gsea_1_cut[,c(2)]
  results_df<-paste0(mir_results_file, '_', Category)
  write.csv(mieaa_all_gsea, paste0(mir_results_file, '_', '.csv' ))
  
  
  ### Convert and return enrich result 
  mieaa_gsea_1$P.adjusted<-as.numeric(mieaa_gsea_1$P.adjusted)
  
  
  mieaa_gsea_1$keyColname=mieaa_gsea_1$Subcategory
  mieaa_gsea_1_ord=mieaa_gsea_1[order(mieaa_gsea_1$P.adjusted),]
  #mieaa_gsea_1_ord_prob<-mieaa_gsea_1_ord
  
  enr_full <- multienrichjam::enrichDF2enrichResult(as.data.frame(mieaa_gsea_1_ord),
                                                    keyColname =  'Subcategory',
                                                    geneColname ='miRNAs.precursors',
                                                    pvalueColname = 'P.adjusted', 
                                                    pvalueCutoff = pvalueCutoff)
  
  
  
  return(list(mieaa_gsea_1, enr_full))
  
}


get_enrich_result_pcgse<-function(all_fs_merged2_pval2){
  #

      all_fs_merged2_pval2$ID<-all_fs_merged2_pval2$Description
        enr_full <- multienrichjam::enrichDF2enrichResult(all_fs_merged2_pval2,
                                                        keyColname =  'Description',
                                                        geneColname ='Description',
                                                        pvalueColname = 'pvalue')
                                                       # pvalueCutoff = pvalueCutoff)
 return(enr_full)
#}
}


get_pval_text<-function(gse, pvalueCutoff_sig){
  text_p1=ifelse(run_ORA,paste0('\n DE: ',  length(gene_list_ora)), '')
  text_p2<-paste0('\n p-adj.< ', pvalueCutoff_sig,': ', length(which(gse@result$p.adjust<pvalueCutoff_sig)), 
                  '\n p-val.< ', pvalueCutoff_sig,': ', length(which(gse@result$pvalue<pvalueCutoff_sig))  )
  text_p=paste0(text_p1, text_p2)
}






###########
get_long_mir_targets<-function(mieaa_targets){
  #' function for mirs 
  #' @param name description
  #'
  #'
  
  dt<-data.table(mieaa_targets)
  
  all_targets_wide<-dcast(dt[, {x1 <- strsplit(miRNAs.precursors, "\\; "); c(list(unlist(x1)), 
                                                                             .SD[rep(seq_len(.N), lengths(x1))])}], Subcategory + miRNAs.precursors ~ V1, length)
  
  all_targets_wide$miRNAs.precursors<-NULL
  all_targets_long<-melt(all_targets_wide)
  all_targets_long_true<-all_targets_long[all_targets_long$value==1, ]
  return(all_targets_long_true)
  
  
}




adjust_unwanted_variance<-function(vsd ){
  #'
  #' Adjust vsd data for unwanted covariates : usable bases and plate 
  #' @param vsd  description
  #' @return vsd_cor corrected vsd
  #'
  colData(vsd)$Usable_Bases_SCALE  <-as.numeric(scale(colData(vsd)$Usable.Bases....))
  colData(vsd)$Plate
  
  retainedCovariates<-colData(vsd)[,c('COHORT')]
  removedCovariates<-colData(vsd)[,c('Usable_Bases_SCALE', 'Plate', '')]
  
  ### Asjustment works on gaussian data so insert vsd or log cplm
  adjusted_data<-empiricalBayesLM(t(as.matrix(assay(vsd))), removedCovariates=removedCovariates, 
                                  retainedCovariates = retainedCovariates )
  
  vsd_cor<-t(adjusted_data$adjustedData)
  return(vsd_cor)
  
}



################ MOFA ####

get_highly_variable_matrix<-function(prefix, VISIT_S, min.count, sel_coh_s,sel_subcoh_s, TOP_N ){
  #''
  #' loads vsd and filters on the fly to save memory! 
  #' @param prefix
  #' 
  #' 
  param_str_tmp<-paste0(prefix, VISIT_S, '_', min.count, '_coh_', sel_coh_s, '_', sel_subcoh_s )

  if (ruv){
    
    # Remove unwanted variance from the vsd data associated with plate and removable bases 
    # 
    # load all and filter 
    # load the corrected dataset - correction is done with all batches together
      print(paste(prefix, ' remove variance'))
      # we do not need the config here but double check that the right file is used.. 

      vst_cor_all_vis_to_load<-paste0(output_files, prefix, 'cell_corr_', cell_corr_mofa, 'vst_cor') # 
      vsd_cor_l=loadRDS(vst_cor_all_vis_to_load) # load the corrected 
      
      ## filter for the specified visit 
      vsd_cor_visit_filt<-filter_se(vsd_cor_l, VISIT = VISIT, sel_coh = sel_coh, sel_sub_coh = sel_subcoh)
      vsd_mat=assay(vsd_cor_visit_filt)





  } else{
  # TODO: Correct all visits together???
  # uncorrected  
  # also defined in config--> check if updated
  deseq_file <-paste0(output_files, param_str_tmp, 'deseq.Rds'); deseq_file
  
  datalist=loadRDS(deseq_file)
  #ddsSE=datalist[[1]]
  vsd=datalist[[2]]
  vsd_mat=assay(vsd)
  dim(vsd)
  }
  # TODO: for mirs load the already normalized and add log2(mirs_expr_norm+1)
  # Perform correction 
  highly_variable_genes_mofa<-selectMostVariable(vsd_mat, TOP_N)
  print(paste('Loaded highly variables files with settings: ', param_str_tmp, TOP_N))
 
 
  
  
  return(highly_variable_genes_mofa)
  
}


prepare_multi_data<-function(p_params, param_str_g_f, param_str_m_f, TOP_GN, TOP_MN, mofa_params){
  #### Takes in the parameters of the input files and loads them 
  #' return: data_full: a list with the 3 modalities 
  # TODO: simplify the reading and setting the feature column to null? 
  #' @param  p_params, param_str_g, param_str_m : these are set by the config.R
  #' 
  
  
  proteins_outfile = paste0(output_files, p_params_out , '_vsn.csv')
  proteins_outfile_csf = paste0(output_files, p_params_csf , '_vsn.csv')
  proteins_outfile_plasma = paste0(output_files, p_params_plasma , '_vsn.csv')
  
  
  proteins_vsn_mat<-as.matrix(read.csv2(proteins_outfile, row.names=1, header=TRUE, check.names = FALSE))
  proteins_csf_vsn_mat<-as.matrix(read.csv2(proteins_outfile_csf, row.names=1, header=TRUE, check.names = FALSE))
  proteins_plasma_vsn_mat<-as.matrix(read.csv2(proteins_outfile_plasma, row.names=1, header=TRUE, check.names = FALSE))


  highly_variable_proteins_mofa<-  selectMostVariable(proteins_vsn_mat, TOP_PN)
  highly_variable_proteins_mofa_plasma<-selectMostVariable(proteins_plasma_vsn_mat, TOP_PN)
  highly_variable_proteins_mofa_csf<-selectMostVariable(proteins_csf_vsn_mat, TOP_PN)
  
  
  ### Start loading mofa data
  proteomics<-as.data.frame(highly_variable_proteins_mofa)
  proteomics_plasma<-as.data.frame(highly_variable_proteins_mofa_plasma)
  proteomics_csf<-as.data.frame(highly_variable_proteins_mofa_csf)
  
  ##### Load mirnas + RNAs 
  ### we use data.table because there are duplicate samples? 
  ### problem with saving of rownmaes 
  highly_variable_mirnas_mofa = get_highly_variable_matrix(prefix='mirnas_', VISIT_S = VISIT_S ,min.count = MIN_COUNT_M, 
                                                        sel_coh_s = sel_coh_s, sel_subcoh_s = sel_subcoh_s, TOP_N=TOP_MN)

  #dim(highly_variable_genes_mofa)
  # EITHER input to vst or put as is normalized
  miRNA<-as.data.frame(highly_variable_mirnas_mofa)
  ##### Load RNA seq: 
  highly_variable_genes_mofa<- get_highly_variable_matrix(prefix='rnas_', VISIT_S = VISIT_S ,min.count = MIN_COUNT_G, 
                             sel_coh_s = sel_coh_s, sel_subcoh_s = sel_subcoh_s, TOP_N=TOP_GN)
  
  #dim(highly_variable_genes_mofa)
  RNA<-as.data.frame(highly_variable_genes_mofa)
  head(rownames(miRNA)); head(colnames(RNA))
  
  data_full<-list(miRNA=as.matrix(miRNA), 
                  RNA=as.matrix(RNA),
                  proteomics=as.matrix(proteomics), 
                  proteomics_plasma=as.matrix(proteomics_plasma), 
                  proteomics_csf=as.matrix(proteomics_csf))
  
  
  return(data_full)
  
  
}




create_multi_experiment<-function(data_full, combined_bl){
  
  ### take a list of three and create the multiassay experiment 
  #' also take metadata and align
  #' @param data_full:  list of threematrices 
  #' @param metadata:
  #' @return mofa_multi, multi assay experiment
  
  RNA= data_full[['RNA']]
  miRNA= data_full[['miRNA']]
  proteomics= data_full[['proteomics']]
  
  proteomics_csf= data_full[['proteomics_csf']]
  proteomics_plasma= data_full[['proteomics_plasma']]
  
  
  assay_full=c(rep('RNA', dim(RNA)[2]),
               rep('miRNA', dim(miRNA)[2]),
               rep('proteomics_csf', dim(proteomics_csf)[2]),
                rep('proteomics_plasma', dim(proteomics_plasma)[2]))

  
  colname = c(colnames(RNA), colnames(miRNA), colnames(proteomics))
  colname = c(colnames(RNA), colnames(miRNA), colnames(proteomics_csf), colnames(proteomics_plasma))
  
  primary=colname
  sample_map=DataFrame(assay=assay_full, primary=primary, colname=colname)
  common_samples_in_assays=unique(colname)
  ### TODO: is it a problem for duplicates when i make PATNO_EVENT_ID the key column? 
  ### Note: HERE WE lost duplicate metadata ie. double clinical measures for one patient
  combined_bl<-combined_bl[!duplicated(combined_bl$PATNO_EVENT_ID),]
  metadata_filt<-combined_bl[match(common_samples_in_assays, combined_bl$PATNO_EVENT_ID),]
  metadata_filt_unique<-metadata_filt[!duplicated(metadata_filt$PATNO_EVENT_ID),]
  metadata_filt_unique<-metadata_filt_unique[!is.na(metadata_filt_unique$PATNO_EVENT_ID),]
  
  metadata_filt_unique$primary<-metadata_filt_unique$PATNO_EVENT_ID
  
  rownames(metadata_filt_unique)<-metadata_filt_unique$PATNO_EVENT_ID
  
  #rownames(metadata_filt)=metadata_filt$PATNO_EVENT_ID
  
  
  mofa_multi<-MultiAssayExperiment(experiments=data_full,
                                   colData = metadata_filt_unique, 
                                   sampleMap=sample_map)
  
  
  return(mofa_multi)
}


filter_se_byExpr<-function(se_filt, min.count=10){
  #'
  #' filter and return without low expression genes 
  #' @param
  #'
  #'
      raw_counts <- assay(se_filt)
      
      idx <- edgeR::filterByExpr(raw_counts,min.count=min.count, group=se_filt$COHORT )
      idx <- edgeR::filterByExpr(raw_counts,min.count=min.count )
      
      raw_counts <- as.matrix(raw_counts[idx, ])
      
      se_filt = se_filt[idx,];
      dim(se_filt)
      return(se_filt)
      
}


create_hist<-function(df, name){
  
  dfm<-melt(df)
  
  p1<-ggplot(dfm, aes(x=value))+ geom_histogram()+ labs(title='mirnas')
  ggsave(paste0(outdir, 'data_histograms',name,  '.jpeg' ), width = 10, height=8)
}



### EVALUATION PURPOSES ######
######## standardize go names 


standardize_go_names<-function(descriptions){
  #'
  #'
  #'
  descriptions=gsub('-', ' ', tolower(descriptions))
  #descriptions=gsub('^[:alnum:]', '', tolower(descriptions))
  
  
  descriptions=gsub("\\'", '', tolower(descriptions))
  descriptions=gsub("\\,", '', tolower(descriptions))
  
  descriptions=gsub('\\(', '', tolower(descriptions))
  descriptions=gsub('\\)', '', tolower(descriptions))
  descriptions=gsub('\\/', '', tolower(descriptions))
  
  descriptions=tolower(descriptions)
  return(descriptions)
  
}







### PREDICTIONS ON SUMMARIZED EXPERIMENT 
## predict on each or on multi assay?

clipping_values<-function(x){
  #'
  #'
  #' @param 
  #'
  higher_val<-quantile(x, 0.99, na.rm=TRUE)
  lower_quant<-quantile(x, 0.02, na.rm=TRUE)
  non_zero_min=min(x[which(x>0)], na.rm = TRUE)
  
  lower_val<-ifelse(non_zero_min>lower_quant,non_zero_min,lower_quant)
  lower_val
  higher_val
  x[which(x<lower_val)]=as.numeric(lower_val)
  
  x[which(x>higher_val)]=as.numeric(higher_val)
  return(x)
  
}



clip_outliers<-function(df1){
  #'
  #' @param 
  #'
  #'
  df1.quantiles <- apply(df1, 1, function(x, prob=0.99) { quantile(x, prob, names=F) })
  
  for (i in 1:dim(df1)[1]){
    df1[i,][ df1[i,]> df1.quantiles[i] ]<- df1.quantiles[i]
  }
  
  return(df1)
}



### time utils 



preprocess_visit<-function(se_filt_V, common,feat_names=NULL, sel_cohorts, clinvars_to_add=c()){
  # 1. Select PD only 
  # 2. Subselected common samples - for training we can use all of them but 
  # for testing only the ones that we have 
  # 3. CPM or VSN
  # 4. Clip outliers 
  # 5. Add patient number and event 

  
  clinvars_to_add<-unique(c('PATNO', 'PATNO_EVENT_ID', 'AGE', 'SEX',  'COHORT',clinvars_to_add))
  
  
  se_filt_V_pd<-se_filt_V[,se_filt_V$COHORT %in% c(sel_cohorts)]
  se_filt_V_pd<-se_filt_V[,se_filt_V$COHORT %in% sel_cohorts]
  se_filt_V_pd<-se_filt_V_pd[,se_filt_V_pd$PATNO %in% common]
  # CPM or VSN? # cpm for plotting, vsn for 
  df_v<-cpm(assay(se_filt_V_pd),  normalized.lib.sizes=TRUE, log=TRUE )
  
  
  # VSN OPTION 
  #df_v<-cpm(assay(se_filt_V_pd),  normalized.lib.sizes=TRUE, log=TRUE )
  
  
  # ??? Why are we clipping? 
  df_v<- clip_outliers(df_v)
  rownames(df_v) = gsub('\\..*', '',rownames(df_v))
  if (is.null(feat_names)){
    df_V_ens=t(df_v)
    
    
  }else{
    df_V_ens<-t(df_v[rownames(df_v) %in% feat_names,])
    
  }
  v_ens=data.frame(df_V_ens)
  v_ens = cbind(v_ens, colData(se_filt_V_pd)[,
      clinvars_to_add])


  
  return(v_ens)
  
}



preprocess_visit_predict<-function(se_filt_V, common,  sel_cohorts, clinvars_to_add=c()){
  # 1. Select PD only 
  # 2. Subselected common samples - for training we can use all of them but 
  # for testing only the ones that we have 
  # 3. CPM or VSN
  # 4. Clip outliers 
  # 5. Add patient number and event 
  #se_filt_V<-se_filt_V[,se_filt_V$COHORT == 1]
  
  clinvars_to_add<-unique(c('PATNO', 'PATNO_EVENT_ID', 'AGE', 'SEX', 'COHORT', clinvars_to_add))
  
  
  se_filt_V_pd<-se_filt_V[,se_filt_V$COHORT %in% sel_cohorts]
  se_filt_V_pd<-se_filt_V_pd[,se_filt_V_pd$PATNO %in% common]
  # CPM or VSN? # cpm for plotting, vsn for 
  df_v<-cpm(assay(se_filt_V_pd),  normalized.lib.sizes=TRUE, log=TRUE )
  
  # trim min counts before size factor estimation? 
  ddsSE <- DESeqDataSet(se_filt_V_pd, 
                        design =as.formula(~COHORT ))
  ddsSE<-estimateSizeFactors(ddsSE)
  
  
  ### separate vsd? 
  # se_filt[]
  # vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
  vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
  
  df_v<-assay(vsd)
    
    
  df_v<- clip_outliers(df_v)
  #df_V_ens<-t(df_v[rownames(df_v) %in% feat_names,])
  v_ens=data.frame(df_v)
  v_ens = cbind(v_ens, colData(se_filt_V_pd)[,clinvars_to_add
            ])
  
  
  return(v_ens)
  
}



create_visits_df<-function(se, clinvars_to_add, feat_names=feat_names, filter_common=TRUE){
  #' create a merged dataframe that includes all visits together 
  #' use a function to help with memory limit 
  #' @param se 
  #' @param
  se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
  se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)
  se_filt_V06<-filter_se(se, VISIT='V06', sel_coh,sel_ps)
  se_filt_V04<-filter_se(se, VISIT='V04', sel_coh,sel_ps)
  
  #Reduce(intersect, list(a,b,c))
  if (filter_common){
    common=intersect(se_filt_V08$PATNO,se_filt_BL$PATNO )
    
  }else{
    common=unique(c(se_filt_V08$PATNO, se_filt_BL$PATNO))
  }
  
  
  ## DO NOT FILTER THE FEAT NAMES HERE 
  v6_ens<-preprocess_visit(se_filt_V06, common=common,feat_names = feat_names,  sel_cohorts = c(1,2), clinvars_to_add =clinvars_to_add )
  v8_ens<-preprocess_visit(se_filt_V08, common=common,feat_names=feat_names, sel_cohorts = c(1,2), clinvars_to_add=clinvars_to_add)
  v4_ens<-preprocess_visit(se_filt_V04, common=common,feat_names=feat_names,  sel_cohorts = c(1,2), clinvars_to_add=clinvars_to_add)
  bl_ens<-preprocess_visit(se_filt_BL, common=common, feat_names=feat_names, sel_cohorts = c(1,2 ), clinvars_to_add=clinvars_to_add)
  
  ######### PLOT molecular markers 
  ### MELT and MERGE 
  v8_melt<-reshape2::melt(v8_ens,id.vars=c(clinvars_to_add) )
  v6_melt<-reshape2::melt(v6_ens,id.vars=c(clinvars_to_add))
  v4_melt<-reshape2::melt(v4_ens,id.vars=c(clinvars_to_add))
  bl_melt<-reshape2::melt(bl_ens,id.vars=c(clinvars_to_add))
  
  
  bl_melt$VISIT<-'BL'
  v4_melt$VISIT<-'V04'
  v8_melt$VISIT<-'V08'
  v6_melt$VISIT<-'V06'
  
  
  
  merged_melt_orig_1<-rbind(bl_melt, v4_melt)
  merged_melt_orig_1<-rbind(merged_melt_orig_1,v6_melt)
  merged_melt_orig_1<-rbind(merged_melt_orig_1,v8_melt)
  
  return(merged_melt_orig_1)
}





###### CLINVARS ####
#df<-df_mofa




get_changes<-function(df,colData_change, t1, t2 ){
  
  #' scale and get change! 
  #'
  
  df_num_1<-as.data.frame(apply(df[paste0(colData_change,'_', t1 )], 2, as.numeric))
  df_num_2<-as.data.frame(apply(df[paste0(colData_change, '_', t2)], 2, as.numeric))
  
  
  colnames(df_num_1)=colData_change
  colnames(df_num_2)=colData_change
  
  df_change<-calc_zscore_change(df_num_1, df_num_2, t2)
 
  
  return(df_change)
}

tBL='BL'
tF1='V13'
tF2='V14'
#t1='BL'
#t2='V13'
#t3='V14'

#colData_change<-scale_vars_diff
#df_all<-fetch_metadata_by_patient_visit(patno_event_ids , combined=combined_bl_log)
#df<-df_all


get_av_change<-function(df,colData_change, tBL,tF1,tF2 ){
  #'
  #' @param  
  #' get average endpoint
  #'
  df_num_BL<-as.data.frame(apply(df[paste0(colData_change,'_', tBL )], 2, as.numeric))
  
  df_num_13<-as.data.frame(apply(df[paste0(colData_change, '_',tF1 )], 2, as.numeric))
  df_num_14<-as.data.frame(apply(df[paste0(colData_change, '_', tF2)], 2, as.numeric))
  
  
  X=list(df_num_13, df_num_14)
  Y <- apply(do.call(cbind, X), 2, as.numeric)
  Y[,c('NP2PTOT_V13','NP2PTOT_V14' ) ]
  
  Y <- array(Y, dim=c(dim(X[[1]]), length(X)))
  av_TF1_TF2=as.data.frame(apply(Y, c(1, 2), mean, na.rm = TRUE))
  colnames(av_TF1_TF2)<-colData_change
  colnames(df_num_BL)<-colData_change
  
  
  
  df_change<-calc_zscore_change(df_num_BL, av_TF1_TF2, paste0(tF1,'_', tF2))
  
  return(df_change)
  
}



calc_zscore_change<-function(df_num_1, df_num_2, t2){
  #'
  #' @param df_change
  #'
  #'
  
  
  
  df_num_bind<-rbind(df_num_1,df_num_2 )
  
  
  ### Get centre and sd from both dataframes - 
  # TODO: FUNCTION
  scaled_attrs1 <- apply(df_num_bind, 2, function(x){
    sx<-scale(x);cn<-attr(sx, 'scaled:center');
    return(cn)
  })
  
  scaled_attrs2 <- apply(df_num_bind, 2, function(x){
    sx<-scale(x)
    sd<-attr(sx, 'scaled:scale')
    return(sd)
  })
  
  df_num_1_scaled <- scale(df_num_1, center=scaled_attrs1, scale=scaled_attrs2)
  df_num_2_scaled <- scale(df_num_2, center=scaled_attrs1, scale=scaled_attrs2)
  df_change=data.frame(df_num_2_scaled-df_num_1_scaled)
  df_change_perc=data.frame(df_num_2_scaled-df_num_1_scaled)/data.frame(df_num_1_scaled+df_num_1_scaled)
  
  
  colnames(df_change) = paste0(colnames(df_change), '_diff_', t2)
  colnames(df_change_perc)=paste0(colnames(df_change_perc), '_diff_', t2, '_perc')
  df_change=cbind(df_change, df_change_perc)
  return(df_change)
}










############## CLUSTERS 


library('dplyr')

get_clinical_clusters_kml<-function(combined_bl_log_sel_pd,y, nbCluster=4, scale_mat=FALSE){
  
  #combined_bl_log_sel_pd=combined_bl_log_sel_pd_to_clust
  #combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel[,'INEXPAGE']=='INEXPD',]
  unique(combined_bl_log_sel_pd$EVENT_ID)
  
  clin_traj<-combined_bl_log_sel_pd[,c('PATNO','EVENT_ID', y)]
  
  clin_traj<-clin_traj[!is.na(clin_traj$EVENT_ID),]
  
  unique(clin_traj$EVENT_ID)
  
  
  clin_traj_wide<-reshape(clin_traj, idvar='PATNO', timevar='EVENT_ID', direction='wide')
  rownames(clin_traj_wide)<-clin_traj_wide$PATNO
  #clinical_clusters<-kmeans((na.omit(clin_traj_wide)[, -1]), centers=centers)
  
  #return(clinical_clusters$cluster)
  
  #install.packages('kml')
  
  ### Clinical trajectory means 
  # REMOVE columns full of NA
  clin_traj_mat<-as.data.frame(sapply((clin_traj_wide)[, -1], as.numeric))
  # ALSO SCALE
  #clin_traj_mat<-clin_traj_mat[!is.na(clin_traj_mat$EVENT_ID),]
  df<-clin_traj_mat
  clin_traj_mat <- as.matrix(df[,colSums(is.na(df))<nrow(df)])
  if (scale_mat){
    clin_traj_mat<-scale(clin_traj_mat)
    
  }
  #devtools::install_github("JimMcL/trajr")
  #library('trajr')
  #trj <- TrajGenerate(200, random = TRUE, angularErrorSd = .25)
  
  #smoothed<-TrajSmoothSG(clin_traj_mat[1,],3,31)
  
  
  
  CLD = kml::cld(clin_traj_mat, timeInData = 1:dim(clin_traj_mat)[2], maxNA = 2)
  #length(CLD)
  #clusters<-kml::kml(CLD, nbRedrawing = 5)
  
  #nbCluster=4
  
  # run choice
  clust_ids<-getClusters(CLD,nbCluster=nbCluster )
  names(clust_ids)<-clin_traj_wide$PATNO
  
  return(clust_ids)
}





get_clinical_clusters<-function(y, centers=4){
  #'
  #' @param 
  #'
  #'
  combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel$COHORT==1,]
  clin_traj<-combined_bl_log_sel[,c('PATNO','EVENT_ID', y)]
  
  clin_traj<-clin_traj[!is.na(clin_traj$EVENT_ID),]
  #clin_traj$months<-unlist(EVENT_MAP[clin_traj$EVENT_ID], use.names = FALSE)
  
  clin_traj_wide<-reshape(clin_traj, idvar='PATNO', timevar='EVENT_ID', direction='wide')
  rownames(clin_traj_wide)<-clin_traj_wide$PATNO
  clinical_clusters<-kmeans((na.omit(clin_traj_wide)[, -1]), centers=centers)
  
  return(clinical_clusters$cluster)
  
}




library('RColorBrewer')
create_venn<-function(venn_list, fname_venn, main){
  
  #######
  #' @param 
  #'
  #'
  myCol2 <- brewer.pal(length(venn_list), "Pastel2")[1:length(venn_list)]
  venn.diagram(venn_list,
               # Circles
               lwd = 2, lty = 'blank', fill = myCol2, cex=2.5,cat.cex=1.5,
               filename =fname_venn, 
               main=main,
               
               output=FALSE)
}






