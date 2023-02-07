
install.packages('minfi')
BiocManager::install('minfiData')

library(data.table)
library(minfi)
library(minfiData)
library(dplyr)
library(stringr)

ppmi_methyl<-read.csv('ppmi/ppmi_data/methylation/ppmi_140_link_list_20210607.csv')
unique_patients<-unique(ppmi_methyl$PATNO)
summary(ppmi_methyl)
## how many instances for each patients??
## how many patients went to BL, V01, V02, V03?
 
colnames(ppmi_methyl)
stats<-table(ppmi_methyl[c('PATNO', 'EVENT_ID' )])
colnames(stats)
stats2<-as.data.frame(stats)
summary(stats)

stats2$Freq<-as.numeric(stats2$EVENT_ID)
ppmi_methyl %>% 
  group_by(EVENT_ID) %>% 
  summarise(n_distinct(PATNO))


ppmi_methyl_120<-read.csv('ppmi/ppmi_data/methylation/beta_post_Funnorm_PPMI_EPICn524final030618_PPMI_120.csv')


### Intersect by visit 

stats_df<-as.data.frame.matrix(stats)

patients_visit_inter<-stats_df[stats_df$BL==1 & stats_df$V04==1  & stats_df$V06==1 & stats_df$V08==1,]
patients_visit_inter<-stats_df[stats_df$BL==1 & stats_df$V04==1& stats_df$V08==1 ,]

row.names(patients_visit_inter)
NROW(patients_visit_inter)

# PATIENTS IDs with methylation data in baseline! 
methyl_BL


# which ones overlap? 
## Proteomics 
files<-list.files(path='ppmi/data/', pattern='Targeted___untargeted_MS-based_proteomics_of_urine_in_PD',
                    full.names = TRUE)
ppmi_prot<-read.csv('Targeted___untargeted_MS-based_proteomics_of_urine_in_PD*.csv')

### Read all proteomics files and bind rowise 
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



### proteomics - TARGETED olink - 4 panels - merge them? 
### Read all proteomics files and bind columnwise 
### bind with more features not needed for now - we are extracting aptients 


files<-list.files(path='ppmi/ppmi_data/proteomics/targeted_olink/plasma/', pattern='*_NPX*',
                  full.names = TRUE)


ppmi_prot<-read.csv(files[[1]])

prot_BL<-ppmi_prot[ppmi_prot$EVENT_ID=='BL',]$PATNO
prot_BL
unique(prot_BL)
stats_df<-as.data.frame.matrix(table(ppmi_prot[c('PATNO', 'EVENT_ID' )]))

patients_visit_inter<-stats_df[stats_df$BL>1 & stats_df$V04>1  & stats_df$V08>1,]
NROW(patients_visit_inter)
rownames(stats_df)

names_split<-sapply(rownames(stats_df),function(x) {
  unlist(strsplit(x,split='\\-'))}
)
patient_number<-sapply(names_split, '[[', 2)
stats_df$PATNO<-patient_number

stats_df_olink<-stats_df




## miRNAs ids 

mirnas<-read.csv2('ppmi/ppmi_data/mirnas/PPMI_sncRNAcounts/all_quantification_matrix_raw.csv/quantification_matrix_raw.final_ids.csv', sep = '\t')

names<-colnames(as.data.frame(mirnas))[-1]

nn<-as.character(names[2])
names_split<-sapply(names,function(x) {
                      unlist(strsplit(x,split='\\.'))}
                                  )

names_split<- strsplit(names,split='\\.')
patient_number<-sapply(names_split, '[[', 3)
visit<-sapply(names_split, '[[', 4)
visit

head(stats_df_olink)
mirnas_patno<-as.data.frame(cbind(patient_number,visit))

stats_df_mirnas<-as.data.frame.matrix(table(mirnas_patno))
stats_df_mirnas$PATNO<-rownames(stats_df_mirnas)


### rnaseq
rnaseq_files<-read.csv('ppmi/ppmi_data/rnaseq/rna_seq_files.txt', row.names = NULL, sep=' ')

stats_df_rnaseq<-as.data.frame.matrix(table(rnaseq_files))
dim(stats_df_rnaseq)
stats_df_rnaseq$PATNO<-rownames(stats_df_rnaseq)




### MERGE



stats_df_rnaseq$PATNO<-rownames(stats_df_rnaseq)

stats_df_olink$PATNO
stats_df_olink<-stats_df_olink %>% rename_at(vars(-PATNO), ~ paste0(., '_p'))
stats_df_mirnas<-stats_df_mirnas %>% rename_at(vars(-PATNO), ~ paste0(., '_m'))
stats_df_rnaseq<-stats_df_rnaseq %>% rename_at(vars(-PATNO), ~ paste0(., '_r'))


colnames(stats_df_rnaseq)
merged<-merge(stats_df_mirnas, stats_df_olink, by = 'PATNO', all=TRUE)
merged_2<-merge(merged, stats_df_rnaseq, by = 'PATNO', all=TRUE)
dim(merged_2)


### RNAs ids
merged_stats<-merged_2


patients_visit_inter<-merged_stats[merged_stats$BL_p>0 & merged_stats$BL_m>0 & 
                                    merged_stats$BL_r>0 &
                                     merged_stats$V04_p>0 & merged_stats$V04_m>0 &
                                     merged_stats$V04_r>0 &
                                     merged_stats$V06_p>0 & merged_stats$V06_m>0 &
                                     merged_stats$V06_r>0,]



patients_visit_inter<-merged_stats[merged_stats$BL_p>0 & merged_stats$BL_m>0 & 
                                     merged_stats$BL_r>0 
                                    ,]



patients_visit_inter<-merged_stats[merged_stats$BL_p>0 & merged_stats$BL_m>0 &
                                     merged_stats$V04_m>0 & merged_stats$V04_p>0 & 
                                     merged_stats$V06_m>0 & merged_stats$V06_p>0 
                                     
                                    ,]

NROW(unique(patients_visit_inter$PATNO))



