
## Prerequisite files  ##

library(dplyr)
library(ggplot2)

source(paste0(script_dir, '/ppmi/extract_metadata.R'))
# deseq2 vst 
#source(paste0(script_dir, 'ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))

library(data.table)
options(rstudio.help.showDataPreview = FALSE) ## RSTUDIO BUG FIX

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)








#### 
combined_choose<-combined[combined$COHORT %in% c(1,2) & combined$INEXPAGE %in% c('INEXPD', 'INEXHC', 'INEXSNCA', 'INEXLRRK2'),]

group_by_var<-'INEXPAGE'
ggplot(combined_choose)+
  geom_density(aes_string(x='AGE_AT_VISIT', fill=group_by_var,
                   color=group_by_var, group=group_by_var), 
               alpha=0.5)



ggplot(combined)+
  geom_histogram(aes(x=AGE_AT_VISIT, fill=as.factor(COHORT), color=COHORT))







graphics.off()

### Create new variables from the averages 
## scopa does not have a total - maybe just add scopa total ? 
### TODO: scopa total needs to be updated : 
# For questions 1-21 (SCAU1 - SCAU21), add 3 points for each response of “9.” 
# For questions 22-25 (SCAU22 - SCAU25), add 0 points for each response of “9.” 

# because the average is the same as np3total
get_totals<-function(combined,sub_pattern, sub_pattern_detect=NULL ){
  #' groups and averages specific columns 
  #' TODO somehwte it is considering NAs as zeros CHECK 
  #' @param sub_pattern
  #' 
  #' 
  #sub_pattern='MCA'
  
  if (is.null(sub_pattern_detect)){
    # definitely works for scau check others 
      #sub_pattern_detect=paste0(sub_pattern,'[1-9]$|', sub_pattern,'[1-9][1-9]$|', sub_pattern,'[1-9][1-9]$')  #For testing
      sub_pattern_detect=paste0(sub_pattern,'[1-9]|', sub_pattern,'[A-Z]')  #For testing
      
    }
  
  
  
    df<-combined[ , grepl( sub_pattern_detect, colnames( combined ) )
                          & !grepl('TOT',  colnames( combined ) ) ]

    df<-as.data.frame(apply(df, 2, as.numeric))
    
    #df
    colnames(df)
    ### If the measure is scopa add 3 if ==9 or 0 for question 22-25
    if (sub_pattern=='SCAU'){
      df=df[, !colnames(df) %in% c('SCAU26A', 'SCAU26AT', 'SCAU26B', 'SCAU26BT', 'SCAU26C', 'SCAU26CT', 'SCAU26D', 'SCAU26DT')]
      cols0<-c('SCAU22','SCAU23','SCAU24','SCAU25')
      df1<-df[, colnames(df) %in% cols0]
      df2<-df[, !(colnames(df) %in% cols0)]
      df1[df1==9]<-0
      df2[df2==9]<-3
      df<-cbind(df2, df1)      

      
      
    } 
    
    
    

    
    # Sum for each patient 
    ind <- rowSums(is.na(df)) == ncol(df)
    df$sca_tot<-rowSums(df, na.rm = TRUE)
    # If all rows are NA then replace zero sum by NA
    #TODO: or do not run in the first place
    df$sca_tot[ind]<-NA 
    
    
    return(df$sca_tot)

}





log_totals<-function(combined, sub_patterns_all){
  #''
  #''
  #''
  
  df<-combined[ , grepl( sub_patterns_all, colnames( combined ) )
                & grepl('_TOT',  colnames( combined ) ) ]
  
  
  df=data.frame(df)
  df=df[sapply(df, is.numeric)]
  df_log=data.frame(sapply(df, function(x) log2(x+1)))
  return(df_log)
  
}

### Up to here is the function to log specific patterns ####


 ppmi_pigd<-function(combined){
   
   tremor<-c('NP2TRMR', 'NP3PTRMR', 'NP3PTRML', 'NP3KTRMR', 'NP3KTRML', 'NP3RTARU', 'NP3RTALU', 'NP3RTARL', 'NP3RTALL', 'NP3RTALJ', 'NP3RTCON')
   pigd<-c('NP2WALK', 'NP2FREZ', 'NP3GAIT', 'NP3FRZGT', 'NP3PSTBL')
   
   tremor_sums=rowMeans(combined[,tremor], na.rm=TRUE)
   pigd_sums<-rowMeans(combined[,pigd], na.rm=TRUE)
   
   TD_PIGD<-tremor_sums/pigd_sums
   TD_PIGD[pigd_sums==0]<-0
   
   length(TD_PIGD)
   combined$TD<-tremor_sums
   combined$PIGD<-pigd_sums
   
   combined$td_pigd_old_on
   
   combined$TD_PIGD_ratio<-TD_PIGD
   if (TD_PIGD>1.15 | TD){
     combined$TD_PIGD<-'TD'
   }
   
   combined$TD_PIGD[combined$TD_PIGD_ratio> 1.15 | (combined$TD>0 & combined$PIGD==0)]<-'TD'
   combined$TD_PIGD[combined$TD_PIGD_ratio<=0.9]<-'PIGD'
   combined$TD_PIGD[combined$TD_PIGD_ratio>0.9 & combined$TD_PIGD_ratio<1.3 ]<-'In'
   
     
     
   combined$TD_PIGD
   table(combined[combined$EVENT_ID=='V06', c('td_pigd_old','TD_PIGD')])
   
  
}

ppmi_rbdsq<-function(combined){
 add1<-c( 'DRMVIVID', 'DRMAGRAC', 'DRMNOCTB',
  'SLPLMBMV', 'SLPINJUR', 'DRMVERBL',
  'DRMFIGHT', 'DRMUMV','DRMOBJFL',
  'MVAWAKEN', 'DRMREMEM', 'SLPDSTRB', 
  'STROKE', 'HETRA', 'PARKISM', 'RLS',
  'NARCLPSY', 'DEPRS', 'EPILEPSY',
  'CNSOTH')
 add1=paste0(add1, '_rbd')
 add1 %in% colnames(combined)
  
  rbdsq<-rowSums(combined[,add1])
  combined$rbdsq<-rbdsq
  curated_total$rem
  combined$rbd
  #combined[, c('RBD_TOT', 'rbdsq', 'rem')]
  
  return(rbdsq)
  
}


combined$rbdsq<-ppmi_rbdsq(combined)
##### Get clinvar changes over time ####




library(stringr)


############



# REM SLEEP BEHAVIOR 
## FIXED 

add_rbd<-c( 'DRMVIVID', 'DRMAGRAC', 'DRMNOCTB',
         'SLPLMBMV', 'SLPINJUR', 'DRMVERBL',
         'DRMFIGHT', 'DRMUMV','DRMOBJFL',
         'MVAWAKEN', 'DRMREMEM', 'SLPDSTRB', 
         'STROKE', 'HETRA', 'PARKISM', 'RLS',
         'NARCLPSY', 'DEPRS', 'EPILEPSY',
         'CNSOTH')
sub_pattern_detect=add_rbd
sub_patterns_rbd<-paste(add_rbd, collapse='|')
rbd_tot<-get_totals(combined=combined, sub_pattern = 'null', sub_pattern_detect = sub_patterns_rbd)

combined$RBD_TOT=rbd_tot
#avs

combined$RBD_TOT



### TODO: STAIAD SHOULD NOT BE ADDED since some of the variables have reverse health outcome 
sub_patterns=c( 'SCAU', 'STAIAD', 'NP3','NP1', 'NP2', 'NP4', 'ESS', 'MCA')
# 1. ADD THE TOTALS TO THE  METADATA  
avs<-sapply(sub_patterns,get_totals, combined=combined)
colnames(avs)<-paste0(colnames(avs), '_TOT')



### TODO: FIX THIS HERE IT CREATES PROBLEM IF  I keep adding and there are duplicates
combined<-cbind(combined,avs)
# add scopa before 
combined$scopa_tot<-combined$SCAU_TOT


## SUBPATTERNS TO LOG 
sub_patterns_2=c(sub_patterns, 'RBD', 'stai_state', 'stai_trait')
sub_patterns_all<-paste(sub_patterns_2, collapse='|')



### from now on work on new vars!! 





#df<-combined[, log_vars]
log_df<-function(df){
  
  df=data.frame(df)
  df=data.frame(apply(df,2, as.numeric))
  df_log=data.frame(sapply(df, function(x) log2(x+1)))
  
  return(df_log)
}

## LOG only the TOTALS !! 
#'
#'
combined$updrs2_score
log_vars<-c('NP3TOT', 'NP2PTOT', 'MCATOT', 'updrs2_score', 'updrs_totscore', 'updrs_totscore_on','updrs2_score', 'updrs3_score', 
            'updrs3_score_on','scopa', 'moca')
df_log2=log_df(combined[, log_vars])
colnames(df_log2)<-paste0(colnames(df_log2),'_LOG')
df_log2$updrs2_score_LOG


df_log<-log_totals(combined,sub_patterns_all = sub_patterns_all)

colnames(df_log)<-paste0(colnames(df_log),'_LOG')
df_log<-cbind(df_log, df_log2)

combined_new<-cbind(combined, df_log)
colnames(df_log)

metadata_output_all<-paste0(output_files, 'combined_log',  '.csv')


combined_new$NP3TOT_LOG

#### ADD FUTURE VISIT #####


### attach future data 



# TODO: also add the outcome total: 
#img_var='NP3_TOT'
## here draw from the original..? 

clinical_scales<-c("NP3TOT" ,  "NP2PTOT"  ,"RBD_TOT",  "MCATOT" ,  "SCAU_TOT", 'NP3TOT_LOG', 'NP2PTOT_LOG', 'updrs3_score', 'updrs2_score', 
                   'updrs3_score_on', 'updrs2_score_LOG', 'NP3TOT_LOG', 'NP2PTOT_LOG', 
                   'updrs3_score_LOG', 'updrs3_score_on_LOG' )
selected_future_vars<-c('PATNO', 'EVENT_ID', 'PDMEDYN', clinical_scales)


cols_fut_visit<-colnames(curated_total_new_cols) # could subselect SOME variables 
patno_event_ids_future<-paste0(combined_new$PATNO, '_', 'V10');
c(cols_fut_visit, selected_future_vars) %in% colnames(combined_new)
combined_future_V10<- fetch_metadata_by_patient_visit(patno_event_ids_future, combined=combined_new)[,c(cols_fut_visit, selected_future_vars)];
'NP3TOT_LOG' %in% colnames(combined_future_V10)

#imaging_variables_diff
patno_event_ids_future<-paste0(combined_new$PATNO, '_', 'V12');
combined_future_V12<- fetch_metadata_by_patient_visit(patno_event_ids_future, combined=combined_new)[,c(cols_fut_visit, selected_future_vars)];


clinical_scales %in% colnames(combined_new)
patno_event_ids_future<-paste0(combined_new$PATNO, '_', 'V14');
combined_future_V14<- fetch_metadata_by_patient_visit(patno_event_ids_future, combined=combined_new)[,selected_future_vars];
patno_event_ids_future<-paste0(combined_new$PATNO, '_', 'V13');
combined_future_V13<- fetch_metadata_by_patient_visit(patno_event_ids_future, combined=combined_new)[,selected_future_vars];




# choose what is available here? 
patno_event_ids_future<-paste0(combined_new$PATNO, '_', 'V16')
combined_future_V16<- fetch_metadata_by_patient_visit(patno_event_ids_future, combined=combined_new)[,selected_future_vars]
combined_future_V16$SCAU_TOT



patno_event_ids_BL<-paste0(combined_new$PATNO, '_', 'BL')
combined_BL<- fetch_metadata_by_patient_visit(patno_event_ids_BL,  combined=combined_new)[,cols_fut_visit]
combined_BL_all<- fetch_metadata_by_patient_visit(patno_event_ids_BL, combined=combined_new)[,c(cols_fut_visit, selected_future_vars)]
combined_BL_all$updrs3_score



# takes a while
# rename and column_bind<

# shall we scale before diff??????? 
# scale by patient though? 
colnames(combined_future_V10)<-paste0(colnames(combined_future_V10), '_V10') # imaging available 
colnames(combined_future_V12)<-paste0(colnames(combined_future_V12), '_V12') # curated available
colnames(combined_future_V16)<-paste0(colnames(combined_future_V16), '_V16')# other variables available
colnames(combined_future_V14)<-paste0(colnames(combined_future_V14), '_V14')# other variables available
colnames(combined_future_V13)<-paste0(colnames(combined_future_V13), '_V13')# other variables available



colnames(combined_BL)<-paste0(colnames(combined_BL), '_BL')# other variables available
colnames(combined_BL_all)<-paste0(colnames(combined_BL_all), '_BL')# other variables available



combined_new<-cbind(combined_new,combined_future_V10 )
combined_new<-cbind(combined_new,combined_future_V12 )
combined_new<-cbind(combined_new,combined_future_V16 )
combined_new<-cbind(combined_new,combined_future_V14 )
combined_new<-cbind(combined_new,combined_future_V13 )

combined_new<-cbind(combined_new,combined_BL_all )

#combined_new$NP3_

write.csv2(combined_new,metadata_output_all, row.names = FALSE)
#combined_bl_log<-combined_new

#combined_new$NP2PTOT_LOG_diff_V14
#combined_new[,'NP2PTOT_LOG_diff_V14']






