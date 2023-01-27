
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
summary(stats)

stats2$Freq<-as.numeric(stats2$EVENT_ID)
ppmi_methyl %>% 
  group_by(EVENT_ID) %>% 
  summarise(n_distinct(PATNO))


# PATIENTS IDs with methylation data in baseline! 
methyl_BL


# which ones overlap? 
## Proteomics
files<-list.files(path='ppmi/data/', pattern='Targeted___untargeted_MS-based_proteomics_of_urine_in_PD',
                    full.names = TRUE)
ppmi_prot<-read.csv('Targeted___untargeted_MS-based_proteomics_of_urine_in_PD*.csv')
read_all<-lapply(files, read.csv)

all_frames<-lapply(read_all, as.data.frame)
all_frames2<-bind_rows(all_frames)

ppmi_prot=all_frames2

prot_BL<-ppmi_prot[ppmi_prot$CLINICAL_EVENT=='BL',]$PATNO

pat_ids<-unique(ppmi_prot$PATNO)

NROW(intersect(prot_BL, methyl_BL))



###INTERSECT for each visit 

events<-unique(ppmi_methyl$EVENT_ID)

sapply(events, function(event_id){
  methyl_visit<-ppmi_methyl[ppmi_methyl$EVENT_ID==event_id,]$PATNO
  prot_visit<-ppmi_prot[ppmi_prot$CLINICAL_EVENT==event_id,]$PATNO
  NROW(intersect(methyl_visit, prot_visit))
})


### RNAs ids


## miRNAs ids 


