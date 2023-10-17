

## which are prodromal to PD
tmp_meta<-samples_metadata(MOFAobject)
tmp_meta$current_cohort<-tmp_meta$COHORT
tmp_meta$DIAG1VIS_num<- as.numeric(gsub('V', '',tmp_meta$DIAG1VIS))
VISIT_num<-as.numeric(gsub('V', '',VISIT))

tmp_meta$DIAG1VIS_num<VISIT_num
# Which patients have already turned to PD? before the specific visit? 
tmp_meta$DIAG1VIS_num<VISIT_num


tmp_meta$current_cohort[(tmp_meta$DIAG1=='PD' & tmp_meta$DIAG1VIS_num<VISIT_num)]<-1
tmp_meta$current_cohort[(tmp_meta$DIAG1=='PD' & tmp_meta$DIAG1VIS_num<VISIT_num)]


tmp_meta$current_cohort

samples_metadata(MOFAobject)$current_cohort<-tmp_meta$current_cohort
