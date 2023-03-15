##### Extract the clinical variables 


# 1. Subject characteristics 
library(data.table)

input_data<-('ppmi/ppmi_data/')
output_files<-'ppmi/output/'

all_files<-list.files(paste0(input_data, 'characteristics/Medical/'), full.names = TRUE)
cl_1<-lapply(all_files, read.csv)

all_files<-list.files(paste0(input_data, 'characteristics/_Subject_Characteristics/'), full.names = TRUE)
cl_2<-lapply(all_files, read.csv)


all_files<-list.files(paste0(input_data, 'Study_Enrollment/'), full.names = TRUE)
cl_3<-lapply(all_files, read.csv)



all_files<-list.files(paste0(input_data, 'ppmi_online/'), full.names = TRUE)
cl_4<-lapply(all_files, read.csv)


lapply(all_files, grep, pattern='')

clinical<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Age_at_visit.csv'))
# todo: add more motor mds-updrs
motor_assess<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_I.csv'))
motor_assess_II<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS_UPDRS_Part_II__Patient_Questionnaire.csv'))
motor_assess_III<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_III.csv'))
motor_assess_IV<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_IV__Motor_Complications.csv'))

non_motor<-read.csv(paste0('ppmi/ppmi_data/SCOPA-AUT.csv'))


VISIT='BL'


### Eventually merge both by visit and patient


motor_assess_all<-merge(motor_assess, motor_assess_III, by=c('PATNO','EVENT_ID'),  suffixes = c("", '_M3'),  all=TRUE); unique(motor_assess_all$PATNO)
motor_assess_all<-merge(motor_assess_all, motor_assess_II,by=c('PATNO','EVENT_ID'),suffixes = c("", '_M2'), all=TRUE) ;unique(motor_assess_all$PATNO)
motor_assess_all<-merge(motor_assess_all, motor_assess_IV,by=c('PATNO','EVENT_ID'),suffixes = c("", '_M4'), all=TRUE) ;unique(motor_assess_all$PATNO)

combined<-merge(motor_assess_all, non_motor, by=c('PATNO','EVENT_ID'), all=TRUE)
combined<-merge(combined, clinical,by=c('PATNO','EVENT_ID'), all=TRUE)

unique(motor_assess_all$PATNO)
#demographics_2<-subset(demographics, select = -c(EVENT_ID))
#combined<-merge(combined, demographics_2,by=c('PATNO'), suffixes = c('.xx', '.de') )


# filtering here only to produce separate files for each visit? 
combined_bl<-combined[combined$EVENT_ID==VISIT,]

## add demographics with suffix if common? 

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
write.csv2(combined,metadata_output_all, row.names = FALSE)


combined

# females would be NA





