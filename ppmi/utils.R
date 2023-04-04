

## Utils 
## Summarized experiment 

### TODO: move to a utils / preprocessing file because it is used also for proteoomics
library(SummarizedExperiment)


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



## Create the summarized experiment by selecting VISITS and cohorts 
filter_se<-function(se, VISIT, sel_coh){
  
  #' Takes the raw file with all counts
  #' Filters summarized experiment by selecting VISITS and cohorts 
  #' @param VISIT
  #' @param sel_coh
  
  ##### 2.   start filtering the experiment  to normalize as appropriate 
  ## Option 1: normalize cohort and EVENT separately!! 
  
  se_filt<-se[,((se$EVENT_ID %in% VISIT) & (se$COHORT %in% sel_coh ))]
  se_filt$EVENT_ID; se_filt$COHORT
  Sample<-colnames(se_filt)
  sample_info<-DataFrame(Sample=Sample)
  
  raw_counts=assays(se_filt)[[1]]
  
  ## filterbyExpr takes cpm so remove from there 
  idx <- edgeR::filterByExpr(raw_counts,min.count=min.count)
  
  length(which(idx))
  raw_counts <- as.matrix(raw_counts[idx, ])
  dim(raw_counts)
  se_filt=se_filt[idx]
  
  ##### Define
  
  ### TODO: Question: Should I input everything into the matrix to normalize? 
  ### And then filter 
  
  ### batch effect and normalization 
  # Create a separate matrix with counts only
  # Include batch information if there is any
  #sample_info$Batch <- as.factor(sample_info$Batch)
  
  
  ### DEFINE THE DESEQ OBJECT with the groups appropriately 
  se_filt$EVENT_ID=as.factor(se_filt$EVENT_ID)
  se_filt$COHORT=as.factor(se_filt$COHORT)
  se_filt$PATNO=as.factor(se_filt$PATNO)
  
  
  return(se_filt)
  
}




### METADTA 

#install.packages('eeptools' )
#library('eeptools')

#as.Date(as.character(new$STATUS_DATE), format = "MM/YY",)

# TODO: fix 
get_age_at_visit<-function(new){
  AGE_AT_VISIT<-as.numeric(gsub('.*/','',new$STATUS_DATE)) - as.numeric(gsub('.*/','',new$BIRTHDT))
  return(AGE_AT_VISIT)
  }
#x_age <- age_calc( as.Date(new$BIRTHDT),          # Convert birth to age
##                   as.Date(new$STATUS_DATE),
#                  units = "years")





#### data specifc 
get_symbols_vector<-function(ens ){
  #' @param ens ensemble ids to conver to symbols 
  #' @returns symbols_ordered the total 
  #'  
  #'  
  
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  symbols <- symbols[!is.na(symbols)]
  symbols_ordered <- symbols[match(ens, names(symbols))]
  na_ind<-is.na(symbols_ordered);
  
  # Add ensembl ids if no symbol found
  symbols_ordered[na_ind]=ens[na_ind]
  return(symbols_ordered)
  
  
}




