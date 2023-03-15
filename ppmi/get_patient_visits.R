
#install.packages('minfi')
#BiocManager::install('minfiData')

library(data.table)
library(minfi)
library(minfiData)

library(dplyr)


library(stringr)
## GMsetEx: GenomicMethylSet, 485,512 features

baseDir <- system.file("extdata", package = "minfiData")
list.files(baseDir)

# tooo sloww
#ppmi_methyl<-read.csv('ppmi/ppmi_data/methylation/ppmi_140_link_list_20210607.csv')
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


#ppmi_methyl_120<-read.csv('ppmi/ppmi_data/methylation/beta_post_Funnorm_PPMI_EPICn524final030618_PPMI_120.csv')


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


prot_files<-list.files(path='ppmi/ppmi_data/proteomics/targeted_olink/plasma/', pattern='*_NPX*',
                  full.names = TRUE)


ppmi_prot<-read.csv(prot_files[[1]])

stats_df<-as.data.frame.matrix(table(ppmi_prot[c('PATNO', 'EVENT_ID' )]))
View(ppmi_prot)

names_split<-sapply(rownames(stats_df),function(x) {
  unlist(strsplit(x,split='\\-'))}
)
patient_number<-sapply(names_split, '[[', 2)
stats_df$PATNO<-patient_number

stats_df_olink_orig<-stats_df

cols<-colnames(stats_df_olink)
library(dplyr)
stats_df_olink_orig<-stats_df_olink_orig %>% 
  mutate(across(!PATNO,  ~1 * (. != 0)))
  


## miRNAs ids 
mirnas_file<-'ppmi/ppmi_data/mirnas/PPMI_sncRNAcounts/all_quantification_matrix_raw.csv/quantification_matrix_raw.final_ids.csv'
mirnas<-read.csv2(mirnas_file, sep = '\t')



names<-colnames(as.data.frame(mirnas))[-1]
names_split<-sapply(names,function(x) {
                      unlist(strsplit(x,split='\\.'))}
                                  )

mirnas_ids<-read.table(mirnas_file, head=TRUE,nrows=1, sep='\t')
names<-colnames(as.data.frame(mirnas_ids))[-1]
names
names_split<- strsplit(names,split='\\.')
patient_number<-sapply(names_split, '[[', 3)
visit<-sapply(names_split, '[[', 4)
visit


mirnas_patno<-as.data.frame(cbind(patient_number,visit))

stats_df_mirnas_orig<-as.data.frame.matrix(table(mirnas_patno))
stats_df_mirnas_orig$PATNO<-rownames(stats_df_mirnas)



### rnaseq
rnaseq_files<-read.table('ppmi/ppmi_data/rnaseq/rna_seq_files.txt', row.names = NULL, sep=' ')
rnaseq_files_IR3<-read.table('ppmi/ppmi_data/rnaseq/IR3_rna_seq_files_ids_sorted.txt', row.names = NULL, sep=' ')
rnaseq_files<-rbind(rnaseq_files,rnaseq_files_IR3 )


stats_df_rnaseq_orig<-as.data.frame.matrix(table(rnaseq_files))
stats_df_rnaseq_orig$PATNO<-rownames(stats_df_rnaseq_orig)
s_df<-stats_df_rnaseq_orig
patients_visit_inter<-s_df[s_df$BL>0 & s_df$V04>0  & s_df$V08>0,]

NROW(patients_visit_inter)
### MERGE

library(dplyr)

stats_df_rnaseq$PATNO<-rownames(stats_df_rnaseq)

stats_df_olink$PATNO
stats_df_olink<-stats_df_olink_orig %>% rename_at(vars(-PATNO), ~ paste0(., '_p'))
stats_df_mirnas<-stats_df_mirnas_orig %>% rename_at(vars(-PATNO), ~ paste0(., '_m'))
stats_df_rnaseq<-stats_df_rnaseq_orig %>% rename_at(vars(-PATNO), ~ paste0(., '_r'))




colnames(stats_df_rnaseq)
merged<-merge(stats_df_mirnas, stats_df_olink, by = 'PATNO', all=TRUE)
merged_2<-merge(merged, stats_df_rnaseq, by = 'PATNO', all=TRUE)
dim(merged_2)


### RNAs ids

# Find groups of patients with intersections 
merged_stats<-merged_2
cols<-c('BL_p', 'BL_r', 'BL_m', 'V04_p', 'V04_r', 'V04_m')
cols<-c('BL_p', 'BL_r', 'BL_m', 'V06_p', 'V06_r', 'V06_m')
cols<-c('BL_p', 'BL_r', 'BL_m', 'V08_p', 'V08_r', 'V08_m')
cols<-c('BL_p', 'BL_r', 'BL_m', 'V04_p', 'V04_r', 'V04_m', 'V06_p', 'V06_r', 'V06_m')
cols<-c('BL_p', 'BL_r', 'BL_m', 'V04_p', 'V04_r', 'V04_m', 'V08_p', 'V08_r', 'V08_m')



subgroup<-merged_stats[, c('PATNO', cols)]
subgroup[subgroup==0]<-NA
NROW(subgroup[complete.cases(subgroup),])

