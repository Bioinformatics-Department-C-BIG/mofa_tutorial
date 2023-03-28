colnames(prot_bl_wide)

combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_', combined$EVENT_ID)
combined[match(colnames(prot_bl_wide),combined$PATNO_EVENT_ID ),]
ids_in_meta<-match(colnames(prot_bl_wide),combined$PATNO_EVENT_ID )
patients_in_meta<-unique(combined$PATNO_EVENT_ID[ids_in_meta])
patients_in_meta<-patients_in_meta[!is.na(patients_in_meta)]
patients_in_meta

# remove duplicates
combined_uniq<-combined[!duplicated(combined$PATNO_EVENT_ID,fromLast=TRUE),]
all(!duplicated(colnames(prot_bl_wide), fromLast=TRUE))





### safe to combine if not duplicates
#### AFTER all the dfs are filtered, order them here
unique(combined$PATNO_EVENT_ID)
combined_uniq$PATNO_EVENT_ID 
ids_ordered<-
combined_match<-combined_uniq[match(patients_in_meta,combined_uniq$PATNO_EVENT_ID ),]
prot_matched<-prot_bl_wide[match(patients_in_meta,colnames(prot_bl_wide))]

dim(cbind(combined_match$PATNO_EVENT_ID, colnames(prot_matched)))


combined_match$COHORT_DEFINITION
