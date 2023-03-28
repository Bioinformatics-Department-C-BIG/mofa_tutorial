

## Utils 
## Summarized experiment 

### TODO: move to a utils / preprocessing file because it is used also for proteoomics
library(SummarizedExperiment)


dim(raw_counts_all)
getSummarizedExperimentFromAllVisits<-function(raw_counts_all, combined){
  #
  raw_counts_all<-raw_counts_all[,!duplicated(colnames(raw_counts_all), fromLast=TRUE)]
  combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)
  
  ### some samples do not exist in metadata so filter them out 
  ## 
  common_samples<-intersect(colnames(raw_counts_all),combined$PATNO_EVENT_ID)
  unique_s<-colnames(raw_counts_all)[!(colnames(raw_counts_all) %in% common_samples)]
  metadata_filt<-combined[match(common_samples, combined$PATNO_EVENT_ID),]
  raw_counts_filt<-raw_counts_all[,match(common_samples, colnames(raw_counts_all))]
  dim(metadata_filt)[1] ==dim(raw_counts_filt)[2]
  
  
  #subset sample names

  se=SummarizedExperiment(raw_counts_filt, colData = metadata_filt)
  
  metadata_filt$COHORT_DEFINITION
  
  return(se)
  
  
}


