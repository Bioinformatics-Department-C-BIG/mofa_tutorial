



library(dplyr)
source(paste0(script_dir, '/ppmi/extract_metadata.R'))
# deseq2 vst 
#source(paste0(script_dir, 'ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))

library(data.table)
options(rstudio.help.showDataPreview = FALSE) ## RSTUDIO BUG FIX

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)


combined$COHORT_DEFINITION

combined[which(combined$NHY==101),]$NHY<-NA






combined$gn

#### 
library(ggplot2)
combined_choose<-combined[combined$COHORT %in% c(1,2,4),]
combined_choose<-combined[combined$COHORT %in% c(1,2) & combined$INEXPAGE %in% c('INEXPD', 'INEXHC', 'INEXSNCA', 'INEXLRRK2'),]

combined$INEXPAGE
group_by_var<-'INEXPAGE'
ggplot(combined_choose)+
  geom_density(aes_string(x='AGE_AT_VISIT', fill=group_by_var,
                   color=group_by_var, group=group_by_var), 
               alpha=0.5)



ggplot(combined)+
  geom_histogram(aes(x=AGE_AT_VISIT, fill=as.factor(COHORT), color=COHORT))




library(ggplot2)


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


df_log<-log_totals(combined,sub_patterns_all = sub_patterns_all)


colnames(df_log)<-paste0(colnames(df_log),'_LOG')
combined_new<-mutate(combined, df_log)
combined_new<-mutate(combined_new, df_log2)

metadata_output_all<-paste0(output_files, 'combined_log',  '.csv')
combined_new$updrs2_score_LOG




#### ADD FUTURE VISIT #####


### attach future data 



# TODO: also add the outcome total: 
#img_var='NP3_TOT'
## here draw from the original..? 

clinical_scales<-c("NP3TOT" ,  "NP2PTOT"  ,"RBD_TOT",  "MCATOT" ,  "SCAU_TOT", 'NP3TOT_LOG', 'NP2PTOT_LOG', 'updrs3_score', 'updrs2_score', 
                   'updrs3_score_on', 'updrs2_score_LOG',
                   'updrs3_score_LOG', 'updrs3_score_on_LOG' )
selected_future_vars<-c('PATNO', 'EVENT_ID', 'PDMEDYN', clinical_scales)


cols_fut_visit<-colnames(curated_total_new_cols) # could subselect SOME variables 
patno_event_ids_future<-paste0(combined_new$PATNO, '_', 'V10');
combined_future_V10<- fetch_metadata_by_patient_visit(patno_event_ids_future, combined=combined_new)[,c(cols_fut_visit, selected_future_vars)];


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


# TODO: make a function to add the scaling by COHORT! 










sel_sam<-MOFAobject@samples_metadata$PATNO_EVENT_ID
sel_pats<-MOFAobject@samples_metadata$PATNO


combined_filt<-combined_new[combined_new$PATNO_EVENT_ID %in% sel_sam,]

# 1. Filter by mofa patients 
combined_filt<-combined_new[combined_new$PATNO %in% sel_pats,]
#combined_filt=combined_new
#### ATTEMPT TO SELECT SAMPLES BECAUSE OF PERCENTILES THAT WE DID NOT USE 
df_log_sel=df_log[combined$PATNO_EVENT_ID %in% sel_sam,]

hist(as.matrix(df_log_sel))

#### Replace variables with the log form 
  
#combined_new<-merge(combined,df_log)
# careful to not fail here
### todo





# just for plotting
combined_new_filt<-combined_new[combined_new$PATNO_EVENT_ID %in% sel_sam,]

combined_new$RBD_TOT

  dim(combined_new)
  
  
  hist(log2(combined_new$stai_state))
  


##################### CLINICAL CHANGES #####################################
  
  #################### 
  # TODO: create function 
#combined %>% 

  #concatenate a patient into one row to add 
  

  library(dplyr)
  library(data.table)


  
  
######################################################



#### Create averages 
##




##### INSPECT DISTIRBUTIONS to decide on appropriate transforms 
###  Create plots only for the specific datasets
### If we are going to normalize/standardize by min-max it is better to do it only for the specific samples? 
i=2
graphics.off()
common_samples=common #### loaded from deseq2_vst_preprocessing script 
combined_p<-combined[combined$PATNO_EVENT_ID %in% common_samples, ]

scales_in_stage<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU_TOT', 'STAIAD_TOT')




####### 







#### Histograms to check the distributions of the clinical variables before and after processing 
for (i in 1:length(scales_in_stage)){
  
    ### 
  
    y=scales_in_stage[i]
    
  
    vals<-combined_p[,y]  
    vals_zero<-vals
    vals_zero[which(vals==0)]<-10^-6
    
    vals_test<-vals[!(is.na(vals) | vals==0)]
    #vals_test<-sample(vals_test, 5000)
    #hist(scale(vals_test, center = FALSE))
    #print(shapiro.test(vals_test))
    
    #hist(combined_to_plot[,y], na.rm=TRUE)
    
    png(paste0(outdir_orig,'metadata/hist_',y,'.png'))
    hist(vals, main=y)
    dev.off()
    
    png(paste0(outdir_orig,'metadata/hist_log_',y,'.png'))
    log_norm<-log2(vals)
    hist(log_norm, main=y)
    dev.off()
    
    log_norm<-log2(vals_zero)
    hist(log_norm)
    
    scaled<-scale(log2(vals_zero), center = FALSE)
    print(summary(scaled))
    hist(scaled)  
    
}





##### AVERAGE ALL THE DATA 
###
###

scaled=DataFrame();scaled_log<-DataFrame();scaled_log_sc<-DataFrame()

for (i in 1:length(scales_in_stage)){
  
  ### scale all the vars to average them 
  ### FOR NP3TP etc. make log 
  
  y=scales_in_stage[i]
  y2<-paste0(y,'_sc')
  comb_scale_data<-combined[,y]
  #print(paste(y, which(comb_scale_data==0)))
  
  comb_scale_data[which(comb_scale_data==0)]<-(10^-6)
  scaled[,y2]<-scale(comb_scale_data, center = FALSE)
  scaled_log[,y2]<-log2(comb_scale_data)
  
  scaled_log_sc[,y2]<-scale(log2(comb_scale_data), center=FALSE)
  
  
  

}

average_if_not_na<-function(df){
  ind <- rowSums(is.na(df)) == ncol(df)
  df$average<-rowMeans(as.data.frame(df), na.rm = TRUE)
  # If all rows are NA then replace zero sum by NA
  #TODO: or do not run in the first place
  df$average[ind]<-NA
  return(df$average)
}


#scaled$average=rowMeans(as.data.frame(scaled), na.rm=TRUE)




























##########################################################
#### skip the above and analyse from file? #####
###### Now filter by relevant patients to make the plots 


### First filter by combined_filt

metadata_output<-paste0(output_files, 'combined_log.csv') 
combined_bl_log<-read.csv2(metadata_output) # combined_bl_log holds the updated data , log, scaled, future visits 


#combined_filt<-combined_new[combined_new$PATNO %in% sel_pats, ]
combined_filt<-combined_bl_log[combined_bl_log$PATNO %in% sel_pats, ]

# inspect patients
#View(combined_filt[combined_filt$PATNO=='3710',])


table(combined_filt$PAG_NAME_M3)

# FILTER OUT non visits and 'R*' - also do this before all the 
combined_filt<-combined_filt[grepl('V',combined_filt$EVENT_ID  ) | grepl('BL',combined_filt$EVENT_ID  ), ]

combined_filt=as.data.frame(combined_filt)


PS_101<-combined[which(combined$NHY==101),]


dim(PS_101[,c('COHORT_DEFINITION','NHY' )])

# REMOVE OUTLIERS FOR plot consistency



combined_filt<-combined_filt %>% 
  dplyr::filter(NHY!=101)
#  filter(PAG_NAME_M3 %in% c('NUPDRS3', 'NUPDRS3A'))

# NUPDRS3A: post dose 






check_dups<-function(combined_to_plot){
  V08_measures<-unique(combined_to_plot[combined_to_plot$EVENT_ID=='V08', c("PATNO","PDSTATE","EVENT_ID", "PD_MED_USE", "NHY_ON", "COHORT")]);
  dup_pats<-V08_measures[(duplicated(V08_measures$PATNO)),]$PATNO
  print(dup_pats)
  return(dup_pats)
}




outl<-max(combined_filt[y], na.rm=TRUE)
group='line_group'
x='EVENT_ID'
shape='PAG_NAME_M3'

#time points
scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU_TOT')


tps<-read.csv(paste0('ppmi/ppmi_data/','visit_tps.csv'), sep=',')
tps#


tps[,1]<-gsub(' ','', tps[,1] )


##  
inds<-match(combined_filt$EVENT_ID, as.character( tps[,1]))
combined_filt$months<-as.numeric(tps[inds,2])
combined$EVENT_ID
x='EVENT_ID'


## Leparates the specific record for each aptient
combined_filt$line_group = with(combined_filt,paste(PATNO,PAG_NAME_M3,PDSTATE, sep='_' ))
combined_filt$line_group = with(combined_filt,paste(PATNO,PAG_NAME_M3,PDSTATE, sep='_' ))



combined_filt$SCAU
combined_to_plot<-combined_filt %>% dplyr::select(c( y, x, group,
                                             scales, 'line_group', 'PATNO' , 'PAG_NAME_M3',
                                'AGE_AT_VISIT', 'COHORT_DEFINITION', 
                                'INEXPAGE', 'PDSTATE', 'PAG_NAME_M4'))


x='months'
combined$INEX
#combined_to_plot$months<-unlist(EVENT_MAP[combined_to_plot$EVENT_ID], use.names = FALSE)

combined_filt$upsit
fw<-'COHORT_DEFINITION'

combined_filt$line_group
### create an average of all clin vars 


group
y='STAGE_LOG_AV'
colour_by<-'PAG_NAME_M4'
colour_by<-'PAG_NAME_M3'
colour_by<-'PDSTATE'

scales<-c('NP1_TOT','NP2_TOT' , 'NP3_TOT', 'NP4_TOT', 'NHY', 'SCAU_TOT', 'PD_MED_USE', 'NP3TOT')
scales<-c( 'NP3_TOT', 'NP2_TOT' )


y='NP1_TOT'
formula_1<-as.formula('~COHORT_DEFINITION')
formula_1<-as.formula('~INEXPAGE')
formula_1<-as.formula('~PDSTATE')


combined_filt$PAG_NAME_M3
combined_to_plot<-combined_filt

combined_to_plot<-combined_to_plot[combined_to_plot$INEXPAGE %in% c('INEXHC', 'INEXPD'),]
combined_to_plot<-combined_to_plot[combined_to_plot$INEXPAGE %in% c( 'INEXPD'),]

combined_to_plot$tau_LLOD


#### sLIDES:  get numbers of unique patients and unique records of medicated/unmedicated at each visit 
sel_pats
SEL_VIS<-'V17'
PATS<-combined_to_plot[combined_to_plot$EVENT_ID==SEL_VIS, c('PATNO', 'PDMEDYN')]
last_visit_patients<-unique(PATS$PATNO)
unique(PATS) %>%
  group_by(PDMEDYN) %>%
  count()



PATS<-combined_to_plot[combined_to_plot$EVENT_ID==SEL_VIS, c('PATNO', 'PDMEDYN', 'PDSTATE')]
unique(PATS) %>%
  group_by(PDSTATE) %>%
  count()


table(as.data.frame(PATS)[, c('PDMEDYN')])




# ensure same samples
## keep pairs - actually could not find any pairs 
combined_to_plot_med<-combined_to_plot[combined_to_plot$PDMEDYN==1,]

pat_dups<-check_dups(combined_to_plot)

V08_measures_paired<-combined_to_plot[combined_to_plot$EVENT_ID=='V08' & combined_to_plot$PATNO %in% pat_dups,]


### IF I FILTER BY PAIRED SAMPLES THEN distributions are NOT normal AND i cannot apply ANOVA 
V08_measures_paired$PDSTATE=as.factor(V08_measures_paired$PDSTATE)
# density plots 
V08_measures_paired_A<-V08_measures_paired[!(V08_measures_paired$PAG_NAME_M3 %in% 'NUPDRS3A'& V08_measures_paired$PDSTATE %in% c('OFF')) , ]

combined_to_plot$PDSTATE


library(dplyr)
df.shapiro_test <- combined_to_plot %>%
  #filter(PAG_NAME_M3%in% 'NUPDRS3' )%>%
  dplyr::filter(PDSTATE%in% c('ON', 'OFF')) %>%
 # filter(PDSTATE %in% 'OFF' )%>%
  mutate(log_np3 = NP3TOT) %>%
  group_by(EVENT_ID) %>%
  nest() %>%
  mutate(t_res =  map(data, ~check_normal(.x) ))
  

combined_to_plot$PDSTATE=as.factor(combined_to_plot$PDSTATE)

df.shapiro_ttest <- combined_to_plot %>%
  #filter(PAG_NAME_M3%in% 'NUPDRS3' )%>%
  filter(PDSTATE%in% c('ON', 'OFF')) %>%
  # filter(PDSTATE %in% 'OFF' )%>%
  mutate(log_np3 = NP3TOT) %>%
  group_by(EVENT_ID) %>%
  filter(n() >= 10)%>%
  nest(-EVENT_ID) %>%
  mutate(log_np3_avg = map(data, ~t.test(NP3TOT~PDSTATE, .x)$p.value))


### RESULTS ON OFF STATE IS IMPORTANT ie. there is a difference ..!!!!
### Check the results and trust only for which they are normal 
t_test_per_visit<-cbind(unlist(df.shapiro_ttest$EVENT_ID),unlist(df.shapiro_ttest$log_np3_avg<0.05))
t_test_per_visit
data=combined_to_plot[combined_to_plot$EVENT_ID%in%'V08',]





  check_normal<-function(data){
    
    
    gg<-data %>%  
      filter(PDSTATE%in% c('ON', 'OFF')) %>%
      
      group_by(PDSTATE) %>%
      filter(n() >= 3)%>%
      mutate(N_Samples = n()) %>%
      nest() %>%
      mutate(Shapiro =  map(data, ~ shapiro.test(.x$NP3TOT)$p.value>0.05))
    
    print(paste('SHAPIRO ',table(data$EVENT_ID), gg$Shapiro))
    # if all normal: 
    if (dim(gg)[1]>2 & all(as.logical(gg$Shapiro))){
      print('check' )
      t_res<-t.test(data$NP3TOT, data$PDSTATE)
      return(t_res)
      
      
    }
  }
  

p<-ggplot(V08_measures_paired,aes(x=log(NP3_TOT), fill=PDSTATE))+
  geom_density(aes(fill=PDSTATE), alpha=0.5)+
  facet_wrap('~PAG_NAME_M3')

p


### ALL THE CONFOUNDERS: PAG_NAME: NUPDRS3 or 3a (POST DOSE)
### PDMEDYN
## Attempt to predict measures from visit and PD state- is it relevant? 
combined_to_plot %>%
  group_by(PATNO, EVENT_ID) %>%
  mutate(n = n()) %>%
  dplyr::filter(n == 2) %>%
  ungroup() %>%
  select(-PATNO, -n)


stats_np3<-combined_to_plot %>% 
  group_by(PATNO, EVENT_ID, PDSTATE ) %>%
  dplyr::filter(EVENT_ID=='V08')%>%
  arrange(PATNO, EVENT_ID)%>%
  summarise(mean=mean(NP3_TOT), sd=sd(NP3_TOT))%>%
  as.data.frame()


stats_np3
ft<-lm(data = combined_to_plot_med, formula = 'NP3_TOT~PDSTATE+EVENT_ID',)
coefficients(summary(ft))


### keep the pairs only 
formula_1<-as.formula('~PAG_NAME_M3')

formula_1<-as.formula('~PDMEDYN')
formula_1<-as.formula('~PD_MED_USE')
formula_1<-as.formula('~PDSTATE')


########## PLOTS FOR SPECIFIC PATIENTS OVER TIME ######################
#######################################################################

### RENAME after we remove the 'R*' VISITS

combined_to_plot$months<-unlist(EVENT_MAP[combined_to_plot$EVENT_ID], use.names = FALSE)

combined_to_plot_med_only<-combined_to_plot[!is.na(combined_to_plot$PD_MED_USE), ]
#### TODO: PLOT SCALES USING SEPARATE GROUPS OF PATIENTS !!! 
## TODO: CGHECK FUTURE 
## TODO: check numbers of patients with molecular data available at future scales 
 scales<-c('NP1_TOT', 'NP3_TOT', 'NP2_TOT', 'NP1_TOT', 'moca')


NROW(unique(combined_to_plot[combined_to_plot$EVENT_ID=='V14','PATNO']))

# TODO: create another filter: plot patients that we have until V14? OR V16 ? TO understand trends 
last_visit_patients 
EVENT_MAP[SEL_VIS]
# 

x='months'
combined_to_plot_last_visit=combined_to_plot[combined_to_plot$PATNO %in% last_visit_patients & (combined_to_plot$months <= EVENT_MAP[SEL_VIS]), ]
#combined_to_plot_last_visit<-combined_to_plot_last_visit[combined_to_plot_last_visit$NP3_TOT_LOG<9,]

combined_to_plot_last_visit$months


combined_to_plot[,c('scopa', 'scopa_tot', 'PATNO', 'PDSTATE', 'REC_ID_SC')]

combined_to_plot_final=combined_to_plot_last_visit
table(combined_to_plot_final$CSFSAA)



### check duplicates 






y='NP3_TOT'
plot_clinvars_by_patient<-function(combined_to_plot_final,x, y, colour_by, shape, facet_var=NULL, line_group='line_group' ){
  #'
  #' @paramcombined_to_plot_final
  #' @x :variable in the x axis--> time
  #' @y : metadata 
  #'
  #'
  #'
  
  # TODO: define the line group here? group by record id? 
  # for scopa tot we do not have 
  
  
  # TODO: pagname m3 does not apply to oteher scales 
  #shape= 
  # fixes error: Error: data must be uniquely named but has duplicate columns
  #5
  
  WIDTH=6
  HEIGHT=3
  
  
 # combined_to_plot_final<-combined_to_plot_final %>% 
#    dplyr::filter(NP3_TOT<500)
  combined_to_plot_final<-combined_to_plot_final%>%
    dplyr::filter(!!as.symbol(y) != '.') %>%
    as.data.frame()

    colnames(combined_to_plot_final) <- make.unique(names(combined_to_plot_final))

  
  
  ### MAKE Numeric to remove dots!
  combined_to_plot_final[, colour_by] = as.factor(combined_to_plot_final[, colour_by])
  combined_to_plot_final[, x]=as.factor(combined_to_plot_final[, x])
  combined_to_plot_final[, y]<-as.numeric(combined_to_plot_final[, y])
  #combined_to_plot_final[, y]<-clip_outliers(combined_to_plot_final[, y], x_times = 5, lower=FALSE)
  combined_to_plot_final[, y]<-clip_outliers(as.data.frame(combined_to_plot_final[, y]))
  
  
  if (y=='urate'){
    combined_to_plot_final<-combined_to_plot_final%>%
      dplyr::filter(!!as.symbol(y) >100) %>%
      as.data.frame()
    
  }
  
  
  p<-ggplot(combined_to_plot_final, aes_string( x=x, color=colour_by, group=line_group))+
    geom_point(aes_string(y=y,color=colour_by, shape=shape))+
    geom_line(aes_string(y=y,color=colour_by, group=line_group),size=0.5, alpha=0.5 )+
    guides( shape='none', group='none')#+
  
  if (!is.null(facet_var) ){
    height=HEIGHT*2
    formula_1=paste0('~', facet_var)
    
    p+facet_wrap(formula_1, nrow = 4)
    p 
    
  }else{
    height=HEIGHT
  }
 
  ggsave(paste0(outdir_orig,'metadata/lines_',SEL_VIS, paste0(facet_var, collapse=''),group,'_', colour_by,'_', y,'.jpeg' ), width=WIDTH, height=height)

  
  p<-ggplot(combined_to_plot_final, aes_string( x=x,y=y))+
    geom_boxplot(aes_string(x=x, fill=colour_by))
  
  
  p

  ggsave(paste0(outdir_orig,'metadata/box_', SEL_VIS, paste0(facet_var, collapse=''),group,'_', colour_by,'_', y,'.jpeg' ), width=WIDTH, height=HEIGHT)
  # ggsave(paste0(outdir_orig,'metadata/violin_', SEL_VIS, paste0(formula_1, collapse=''), y,'.jpeg' ), width=10, height=7)
}



scales<-c('NP1_TOT', 'NP3_TOT', 'NP2_TOT', 'NP1_TOT', 'NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'NP2_TOT_LOG', 'NP3_TOT_LOG')

# , 'NHY', 'td_pigd' --> get counts for these ones at every time point eg. counts of 2, counts of 3 and plot those
colour_by='PD_MED_USE'
facet_var='PDSTATE'
colour_by='PDSTATE'


for (y in scales){
  plot_clinvars_by_patient(combined_to_plot_final,x, y, colour_by=colour_by, shape=shape, facet_var = facet_var  )

}

# scales not taken at off and on
combined_to_plot_final
colour_by='PDMEDYN'

combined_to_plot_final$upsit_pctl
facet_var='COHORT'
combined_to_plot_final[, c('PATNO','PDMEDYN', 'MCA_TOT', 'upsit_pctl')]
scales<-c( 'SCAU_TOT', 'scopa_tot', 'MCA_TOT', 'bjlot', 'moca', 'upsit', 'MSEADLG', 'ESS_TOT')
scales<-c( 'stai' , 'stai_state', 'stai_trait', 'VLTANIM', 'gds', 'hvlt_immediaterecall', 'ess')

for (y in scales){
  plot_clinvars_by_patient(combined_to_plot_final,x, y, colour_by, shape, facet_var=NULL , line_group='PATNO')
  }

graphics.off()
combined_to_plot_final$ptau




combined_to_plot_final[,c('asyn', 'PATNO_EVENT_ID')]
combined_to_plot_final$asyn<-as.numeric(combined_to_plot_final$asyn)
  combined_to_plot_final$asyn
shape
#iMAGING
combined_to_plot_final$CSFSAA

# ADD upsit, semantic fluency 
 
 # un==imaging

scales=c('ptau', 'asyn', 'tau', 'abeta', 'tau', 'CSFSAA', 'nfl_csf', 'total_di_18_1_BMP', 'hemohi', 'urate')
scales=c('DATSCAN_PUTAMEN_L', 'DATSCAN_CAUDATE_L', 'DATSCAN_CAUDATE_R', 'DATSCAN_PUTAMEN_R', 'ips_caudate', 'mean_caudate', 'mean_striatum', 
         'con_caudate', 'con_putamen', 'con_striatum', 'lowput_ratio')
combined_to_plot_final$DBSYN
facet_var=NULL
#facet_var='DBSYN'

combined_to_plot_final[, scales]
for (y in scales){
  plot_clinvars_by_patient(combined_to_plot_final,x, y, colour_by, shape, facet_var = facet_var, line_group='PATNO' )
}

graphics.off()

combined_to_plot_final[combined_to_plot_final$EVENT_ID=='V16', 'mean_striatum']

combined_to_plot_final
##

curated_v1[curated_v1$EVENT_ID=='V12', 'mean_striatum']
curated_v2[curated_v2$EVENT_ID=='V16', 'mean_striatum']
as.numeric(curated_v2[curated_v2$EVENT_ID=='V08', 'urate'])
as.numeric(curated_v1[curated_v1$EVENT_ID=='V08', 'urate'])

as.numeric(combined_to_plot_final[combined_to_plot_final$EVENT_ID=='V08', 'urate'])

curated_v1$urate
curated_v2$mean_striatum

