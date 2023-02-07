##### Extract the clinical variables 


# 1. Subject characteristics 

input_data<-('ppmi/ppmi_data/')
clinical<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Age_at_visit.csv'))
demographics<-read.csv(paste0(input_data, 'characteristics/_Subject_Characteristics/Demographics.csv'))

motor_assess<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_I.csv'))
motor_assess<-read.csv(paste0(input_data, 'motor_assess/Motor___MDS-UPDRS/MDS-UPDRS_Part_I.csv'))
non_motor<-read.csv(paste0('ppmi/ppmi_data/SCOPA-AUT.csv'))


clinical_BL<-clinical[clinical$EVENT_ID=='BL',]
length(unique(clinical_BL$PATNO))


motor_assess_BL<-motor_assess[motor_assess$EVENT_ID=='BL',]
non_motor_BL<-non_motor[non_motor$EVENT_ID=='BL',]

non_motor

# 2. 
