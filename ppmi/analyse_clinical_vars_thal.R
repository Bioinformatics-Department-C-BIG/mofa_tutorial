


### CLINICAL DATA 

grep('Number',colnames(samplesheet) )


library(dplyr)

process_clinvars<-function(samplesheet){
samplesheet <- samplesheet %>%
  mutate(across(starts_with("Total|Number"),
                ~ as.numeric(as.character(.))))
}


