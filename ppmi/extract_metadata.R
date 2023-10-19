##### Extract the clinical variables 


# 1. Subject characteristics 
library(data.table)

input_data<-('ppmi/ppmi_data/')
output_files<-'ppmi/output/'
source('ppmi/utils.R')
#all_files<-list.files(paste0(input_data, 'characteristics/Medical/'), full.names = TRUE)

#cl_1<-lapply(all_files, read.csv)
# only add the ones with event id 

#new_cl1<-cl_1[grepl('EVENT', cl_1 )]


#cl1_merged<-Reduce(function(x, y) merge(x, y, all=TRUE, by=c('PATNO','EVENT_ID')), new_cl1)  


#all_files<-list.files(paste0(input_data, 'characteristics/_Subject_Characteristics/'), full.names = TRUE)
#subject_characteristics<-lapply(all_files, read.csv)


#all_files<-list.files(paste0(input_data, 'Study_Enrollment/'), full.names = TRUE)
#cl_3<-lapply(all_files, read.csv)


#all_files<-list.files(paste0(input_data, 'Study_Enrollment/'), full.names = TRUE)
#cl_3<-lapply(all_files, read.csv)



#all_files<-list.files(paste0(input_data, 'ppmi_online/'), full.names = TRUE)
#cl_4<-lapply(all_files, read.csv)


# these do not have event id
# TODO: genetic consensus
ps<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Participant_Status.csv'))
genetics<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/iu_genetic_consensus_20220310.csv'))
clinical<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Age_at_visit.csv'))
demographics<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Demographics.csv'))


dim(characteristics)
characteristics<-merge(clinical, ps, by=c('PATNO'), suffixes = c("", '_ps'),  all=TRUE)
dim(characteristics)


othfeatpd<-read.csv(paste0(input_data, 'characteristics/Medical/Other_Clinical_Features.csv'))
prodromal_history<-read.csv(paste0(input_data, 'characteristics/Medical/Prodromal_History.csv'))



# todo: add more motor mds-updrs
motor_assess<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_I.csv'))
motor_assess_II<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS_UPDRS_Part_II__Patient_Questionnaire.csv'))

#motor_assess_III<-read.csv2(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_III.csv'), sep=',', stringsAsFactors = FALSE)
motor_assess_III<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_III-short-date.csv'))

motor_assess_IV<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_IV__Motor_Complications.csv'))
#motor_assess_IV<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/'))

motor_assess_III$ONEXAMDT

## cognition####

non_motor<-read.csv(paste0('ppmi/ppmi_data/SCOPA-AUT.csv'))
non_motor_moca<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/Montreal_Cognitive_Assessment__MoCA_.csv'))
non_motor_stait<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/State-Trait_Anxiety_Inventory.csv'))
pr_outcome<-read.csv(paste0('ppmi/ppmi_data/patient_outcomes_predictions.csv'))



### Sleepiness scale 

olfactory<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/University_of_Pennsylvania_Smell_Identification_Test__UPSIT_.csv'))
epworth<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/Epworth_Sleepiness_Scale.csv'))

rbd<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/REM_Sleep_Behavior_Disorder_Questionnaire.csv'))
dat<-read.csv(paste0('ppmi/ppmi_data/imaging/DaTSCAN/DaTScan_Analysis_23Jun2023.csv'))
bent<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/Benton_Judgement_of_Line_Orientation.csv'))
ger<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/Geriatric_Depression_Scale__Short_Version_.csv'))
sem<-read.csv(paste0('ppmi/ppmi_data/Non-motor_Assessments/Modified_Semantic_Fluency.csv'))


### add biospecimen 



curated_v1<-read.csv(paste0('ppmi/ppmi_data/curated data/Curated_Data_Cuts/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv'))
curated_v2<-read.csv(paste0('ppmi/ppmi_data/curated data/PPMI_Curated_Data_Cut_Public_20230612.csv'))
curated_total<-merge(curated_v1, curated, by='PATNO', )
## first bind the common columns 
common_cols<-intersect(colnames(curated_v1), colnames(curated_v2))
curated_total<-rbind(curated_v1[, common_cols],curated[, common_cols])
curated_v2_unique<-curated_v2[, (!(colnames(curated_v2) %in% common_cols))| colnames(curated_v2) %in%c('PATNO', 'EVENT_ID')  ]; 
curated_v1_unique<-curated_v1[, (!(colnames(curated_v1) %in% common_cols))| colnames(curated_v1) %in%c('PATNO', 'EVENT_ID')  ]; 

curated_total_new_cols<-merge(curated_total, curated_v2_unique, by=c('PATNO', 'EVENT_ID'))
curated_total_new_cols<-merge(curated_total_new_cols, curated_v1_unique, by=c('PATNO', 'EVENT_ID'))

colnames(curated_total_new_cols)




### Eventually merge both by visit and patient
## These measures are per patient and event id

motor_assess_all<-merge(motor_assess, motor_assess_III, by=c('PATNO','EVENT_ID'),  suffixes = c("_M1", '_M3'),  all=TRUE); unique(motor_assess_all$PATNO)
motor_assess_all<-merge(motor_assess_all, motor_assess_II,by=c('PATNO','EVENT_ID'),suffixes = c("", '_M2'), all=TRUE) ;unique(motor_assess_all$PATNO)
motor_assess_all<-merge(motor_assess_all, motor_assess_IV,by=c('PATNO','EVENT_ID'),suffixes = c("", '_M4'), all=TRUE) ;unique(motor_assess_all$PATNO)

combined<-merge(motor_assess_all, non_motor, by=c('PATNO','EVENT_ID'), suffixes = c("", '_SC'),  all=TRUE)
combined<-merge(combined, clinical,by=c('PATNO','EVENT_ID'),suffixes = c("", 'cl'), all=TRUE)

### everything was collected at screening phase only ####
combined<-merge(combined, demographics,by=c('PATNO'),suffixes = c("", 'd'), all=TRUE)
combined<-merge(combined, genetics,by=c('PATNO'),suffixes = c("", '_gen'), all=TRUE)



### cognition
combined<-merge(combined, othfeatpd,by=c('PATNO','EVENT_ID'), suffixes = c("", '_oth'),  all=TRUE)
combined<-merge(combined, prodromal_history,by=c('PATNO','EVENT_ID'), suffixes = c("", '_prod'),  all=TRUE)

combined<-merge(combined, non_motor_moca,by=c('PATNO','EVENT_ID'), suffixes = c("", '_moca'),  all=TRUE)

combined<-merge(combined, non_motor_stait,by=c('PATNO','EVENT_ID'), suffixes = c("", '_st'),  all=TRUE)

combined<-merge(combined, epworth,by=c('PATNO','EVENT_ID'), suffixes = c("", '_ep'),  all=TRUE)
rbd_vals<-colnames(rbd)[!colnames(rbd) %in% c('PATNO','EVENT_ID')]
colnames(rbd)[!colnames(rbd) %in% c('PATNO','EVENT_ID')]<-paste0(rbd_vals, '_rbd')



combined<-merge(combined, curated_total_new_cols,by=c('PATNO','EVENT_ID'), suffixes = c("", '_cur'),  all=TRUE)

combined<-merge(combined, rbd,by=c('PATNO','EVENT_ID'), suffixes = c("", '_rbd'),  all=TRUE)
combined<-merge(combined, dat,by=c('PATNO','EVENT_ID'), suffixes = c("", '_dat'),  all=TRUE)
combined<-merge(combined, olfactory,by=c('PATNO','EVENT_ID'), suffixes = c("", '_olf'),  all=TRUE)
combined<-merge(combined, bent,by=c('PATNO','EVENT_ID'), suffixes = c("", '_bn'),  all=TRUE)
combined<-merge(combined, ger,by=c('PATNO','EVENT_ID'), suffixes = c("", '_gr'),  all=TRUE)
combined<-merge(combined, sem,by=c('PATNO','EVENT_ID'), suffixes = c("", '_sem'),  all=TRUE)



## these measures are per patient 
combined<-merge(combined, ps,by=c('PATNO'), suffixes = c("", '_ps'),  all=TRUE)
combined<-merge(combined, pr_outcome,by=c('PATNO'), suffixes = c("", '_pr'),  all=TRUE)

V10_mean_striatum<-curated_total_new_cols[curated_total_new_cols$EVENT_ID=='V10', c('mean_striatum', 'PATNO')]
combined<-merge(combined, V10_mean_striatum,by=c('PATNO'), suffixes = c("", '_V10'),  all=TRUE)

#which(!is.na(combined$Outcome))


#### FIX age and sex
### OUTPUT THE FILTERED se_filt 

ind<-which(is.na(combined$AGE_AT_VISIT))
combined$AGE<-combined$AGE_AT_VISIT
combined[ind,'AGE' ]<-get_age_at_visit(combined[ind,])
## Turn to factors for deseq
combined$SEX<-as.factor(combined$SEX)
combined$AGE_SCALED<-scale(combined$AGE)


combined$OFFPDMEDDT
combined$INFODT_M1
combined$OFFEXAMDT
combined$HRPOSTMED ### hours since last dose 

combined[,c('ONEXAM', 'OFFEXAM','PDMEDYN','ORIG_ENTRY_M3',  'INFODT_M3','NTEXAMDT',  'OFFEXAMDT' ,'OFFEXAMTM', 'OFFPDMEDDT', 'OFFPDMEDTM')]
combined$PATNO
rbd
#demographics_2<-subset(demographics, select = -c(EVENT_ID))
#combined<-merge(combined, demographics_2,by=c('PATNO'), suffixes = c('.xx', '.de') )


combined$SCAU26CT<-tolower(combined$SCAU26CT)
combined$SCAU26CT<-as.factor(combined$SCAU26CT)

levels(as.factor(combined$SCAU26CT))
# filtering here only to produce separate files for each visit? 

## add demographics with suffix if common? 

### Add new features here

combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)
combined$AGE

metadata_output_all<-paste0(output_files, 'combined',  '.csv')


### TODO: after this run the logs in analyse clinical vars 
write.csv2(combined,metadata_output_all, row.names = FALSE)

combined$COHORT_DEFINITION
#View(combined[combined$PATNO=='4125',])
combined$NHY


# females would be NA

MOFAobject@samples_metadata$PATNO

demographics[which(demographics$PATNO==3386),]

demographics$PATNO
common

combined$GBA_POS

## add genetics





