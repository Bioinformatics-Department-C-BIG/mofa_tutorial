


### GET FUTURE SCALES 
input_df=df_mofa
## TODO: create a function: fetch future data using diffs 

df_all<-fetch_metadata_by_patient_visit(input_df$PATNO_EVENT_ID , combined=combined_bl_log)
hist(df_all$NP2PTOT_V16)
hist(df_all$NP2_TOT_V16)
colData_change<-c('updrs3_score', 'con_putamen', 'hi_putamen', 'mean_striatum','updrs2_score', 'moca')
t1<-'BL';  t2='V10';

#df=df_all
df_change1= get_changes(input_df,colData_change, t1, t2 )
colnames(df_change)
df_all<-cbind(df_all, df_change1)
t1<-'BL';  t2='V16'; 
df_all$MCA_TOT_V16
df_all$moca_BL
df_all$moca_V10
df_all$MCA_TOT_V16
colData_change=c('NP3TOT', 'NP2PTOT', 'RBD_TOT', 'MCA_TOT' )
#df_all$NP3_TOT_V16
df_change= get_changes(df_all,colData_change, t1, t2 )
df_all<-cbind(df_all, df_change)
df_all$updrs2_score_diff_V16

df_all<-fetch_metadata_by_patient_visit(input_df$PATNO_EVENT_ID , combined=combined_bl_log, max_m='NP3TOT')
df_change


################# ##############


df_all<-fetch_metadata_by_patient_visit(vsd$PATNO_EVENT_ID , combined=combined_bl_log)
df_all$NP3_TOT_V16

colData_change<-c('updrs3_score', 'con_putamen', 'hi_putamen', 'updrs2_score', 'moca')
t1<-'BL';  t2='V10';

df_change1= get_changes(df_all,colData_change, t1, t2 )
colnames(df_change)
df_all<-cbind(df_all, df_change1)
df_change1<-cbind(df_all$PATNO, df_change1)


t1<-'BL';  t2='V16';
colData_change=c('NP3TOT', 'NP2PTOT', 'RBD_TOT' , 'MCA_TOT')
#df_all$NP3_TOT_V16
df_change2= get_changes(df_all,colData_change, t1, t2 )
df_all<-cbind(df_all, df_change2)
df_change=cbind(df_change1, df_change2)
df_change$PATNO<-df_change$`df_all$PATNO`


### NOEW MERGE WITH MOFA 



############# AFTER WE ADDE THIS TO MOFA
all_diff<-all_diff_variables[all_diff_variables %in% colnames(cors_all)]
cors_all[, all_diff]


