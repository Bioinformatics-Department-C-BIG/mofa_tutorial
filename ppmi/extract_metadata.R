##### Extract the clinical variables 


# 1. Subject characteristics 

input_data<-('ppmi/ppmi_data/')
clinical<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Age_at_visit.csv'))


# todo: add more motor mds-updrs
motor_assess<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_I.csv'))
motor_assess_III<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_III.csv'))

non_motor<-read.csv(paste0('ppmi/ppmi_data/SCOPA-AUT.csv'))



### Eventually merge both by visit and patient
clinical_BL<-clinical[clinical$EVENT_ID=='BL',]
length(unique(clinical_BL$PATNO))


motor_assess_BL<-motor_assess[motor_assess$EVENT_ID=='BL',]
motor_assess_III_BL<-motor_assess_III[motor_assess_III$EVENT_ID=='BL',]

NROW(motor_assess_BL$PATNO)

NROW(unique(motor_assess_BL$PATNO))


non_motor_BL<-non_motor[non_motor$EVENT_ID=='BL',]

motor_assess_BL<-merge(motor_assess_BL, motor_assess_III_BL, by='PATNO')

motor_assess %>% 
    summarise(n_distinct(unlist(across(PATNO:EVENT_ID))))

motor_assess_III_dt$NHY

motor_assess_III_dt<-as.data.table(motor_assess_III)
NROW(motor_assess_III_dt[, .N, by = c('PATNO', 'EVENT_ID')])

NROW(motor_assess)

motor_assess_all<-merge(motor_assess, motor_assess_III, by=c('PATNO','EVENT_ID'))





motor_assess_all<-as.data.table(motor_assess_all)
NROW(motor_assess_all[, .N, by = c('PATNO', 'EVENT_ID')])


combined<-merge(motor_assess_all, non_motor, by=c('PATNO','EVENT_ID'))

combined<-merge(combined, clinical,by=c('PATNO','EVENT_ID'))
demographics_2<-subset(demographics, select = -c(EVENT_ID))
combined<-merge(combined, demographics_2,by=c('PATNO'), suffixes = c('.xx', '.de') )


combined_bl<-combined[combined$EVENT_ID=='BL',]

## add demographics with suffix if common? 



write.csv2(combined_bl,'combined_bl.csv', row.names = FALSE)

# females would be NA
combined_bl$SCAU23

combined_bl$SEX


