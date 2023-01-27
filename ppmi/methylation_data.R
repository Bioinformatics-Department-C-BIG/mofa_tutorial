
install.packages('minfi')
BiocManager::install('minfiData')

library(minfi)
library(minfiData)
## RGsetEx: RGChannelSet, 622,399 features
MsetEx <- preprocessRaw(RGsetEx)
## MsetEx: MethylSet, 485,512 features
GMsetEx <- mapToGenome(MsetEx)
## GMsetEx: GenomicMethylSet, 485,512 features

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)

library(dplyr)
ppmi_methyl<-read.csv('ppmi/data/ppmi_140_link_list_20210607.csv')
unique_patients<-unique(ppmi_methyl$PATNO)
summary(ppmi_methyl)
## how many instances for each aptients??
 
colnames(ppmi_methyl)
stats<-table(ppmi_methyl[c('PATNO', 'EVENT_ID' )])
colnames(stats)
stats2<-as.data.frame(stats)
table(stats2$)
summary(stats)

stats2$Freq<-as.numeric(stats2$EVENT_ID)
ppmi_methyl %>% 
  group_by(EVENT_ID) %>% 
  summarise(n_distinct(PATNO))

methyl_BL<-ppmi_methyl[ppmi_methyl$EVENT_ID=='BL',]$PATNO
# PATIENTS IDs with methylation data in baseline! 
methyl_BL


# which ones overlap? 
## Proteomics




### RNAs ids


## miRNAs ids 


