

library(dplyr)
library(data.table)
library(stringr)


#### Read in the mirnas 
mirnas_rpmmm<-read.csv2('ppmi/ppmi_data/mirnas/PPMI_sncRNAcounts/mirna_quantification_matrix_rpmmm_norm.csv/std_quantification_rpmmm_norm_mirna.final_ids.csv', sep = '\t')
mirnas_rpmmm<-read.csv2('ppmi/ppmi_data/mirnas/PPMI_sncRNAcounts/mirna_quantification_matrix_raw.csv/std_quantification_raw_mirna.final_ids.csv', sep = '\t')

names<-colnames(mirnas_rpmmm)[-1]

#names_split<- strsplit(names,split='\\.')
#visit<-sapply(names_split, '[[', 4)
#bl_visits<-visit=='BL'


mirnas_BL<-select(mirnas_rpmmm,contains("BL"))
mirnas_BL


names_split<- strsplit(names(mirnas_BL),split='\\.')

PATNO<-sapply(names_split, '[[', 3)
length(PATNO)
colnames(mirnas_BL)<-PATNO
rownames(mirnas_BL)<-mirnas_rpmmm$miRNA

length((mirnas_BL))
tail(colnames(mirnas_BL), 400)
rownames(mirnas_BL)



library(tidyr)
##### PROTEOMICS - OLINK
csf=1
output_files<-'ppmi/output/csf/'
if (csf){
  prot_files<-list.files(path='ppmi/ppmi_data/proteomics/targeted_olink/csf/', pattern='*_NPX*',
                         full.names = TRUE)
  outname<-'proteomics_csf_bl.csv'
  outname2<-'proteomics_csf_bl_no_log.csv'
  
}else {
  
  prot_files<-list.files(path='ppmi/ppmi_data/proteomics/targeted_olink/plasma/', pattern='*_NPX*',
                         full.names = TRUE)
  outname<-'proteomics_bl.csv'
  outname2<-'proteomics_bl_no_log.csv'
}


prot_files



# Merge the 4 datasets together 
read_all<-lapply(prot_files, read.csv)
all_frames<-lapply(read_all, as.data.frame)
all_frames2<-bind_rows(all_frames)
ppmi_prot=all_frames2

### remove PPMI- suffix from PATNO column 
ppmi_prot<-ppmi_prot %>% mutate(across('PATNO', str_replace, 
                            'PPMI-', ''
                            ))
## Filter baseline
prot_bl<-ppmi_prot[ppmi_prot$EVENT_ID=='BL',]

## Remove unecessary columns 
prot_bl_matrix<-as.data.frame(subset(prot_bl,
                    select=-c(LOD, update_stamp, OLINKID, UNIPROT, MISSINGFREQ,
                              PANEL, EVENT_ID, PLATEID, INDEX, PANEL_LOT_NR)))
## Extract QC pass only 
prot_bl_matrix<-prot_bl_matrix[prot_bl_matrix$QC_WARNING=='PASS',]

hist(prot_bl_matrix$NPX)



## Reshape: Make a wide matrix with patient in columns 
prot_bl_tbl<-as.data.table(prot_bl_matrix)
prot_bl_wide<-data.table::dcast(prot_bl_tbl,  ASSAY ~ PATNO,
                         value.var ='NPX', fun.aggregate = mean)


prot_bl_wide<-as.data.frame(prot_bl_wide)
rownames(prot_bl_wide)<-prot_bl_wide$ASSAY
prot_bl_wide$ASSAY<-NULL
row.names(prot_bl_wide)


## Remove the log normalization
prot_bl_wide_unlog<-as.data.frame(sapply(prot_bl_wide, function(x){2**(as.numeric(x))}),
                                    row.names =row.names(prot_bl_wide) )


## Write output both log and not logged
fwrite(prot_bl_wide, paste0(output_files, outname), row.names = TRUE)
fwrite(prot_bl_wide_unlog, paste0(output_files, outname2), row.names = TRUE)


df<-prot_bl_wide_unlog[prot_bl_wide_unlog<10]
hist(as.numeric(as.matrix(df)))




##### PROTEOMICS - URINE MS/MS
prot_files<-list.files(path='ppmi/ppmi_data/proteomics/untargeted/', pattern='*',
                       full.names = TRUE)
outname<-'proteomics_untargeted.csv'


# Merge the 4 datasets together 
read_all<-lapply(prot_files, read.csv)
all_frames<-lapply(read_all, as.data.frame)
all_frames2<-bind_rows(all_frames)
ppmi_prot=all_frames2

### remove PPMI- suffix from PATNO column 
ppmi_prot<-ppmi_prot %>% mutate(across('PATNO', str_replace, 
                                       'PPMI-', ''
))
## Filter baseline
prot_bl<-ppmi_prot[ppmi_prot$CLINICAL_EVENT=='BL',]

## Remove unecessary columns 
remove<-c('COHORT', 'TYPE', 'UNITS', 'RUNDATE', 'PROJECTID', 'PI_NAME', 'PI_INSTITUTION', 'update_stamp', 'SEX', 'EVENT_ID')
prot_bl_matrix<-as.data.frame(subset(prot_bl,
                                     select=-c(COHORT, TYPE, UNITS, RUNDATE, PROJECTID, PI_NAME, PI_INSTITUTION, update_stamp, SEX, CLINICAL_EVENT)))


## Extract QC pass only 
#prot_bl_matrix<-prot_bl_matrix[prot_bl_matrix$QC_WARNING=='PASS',]

hist(prot_bl_matrix$TESTVALUE)



## Reshape: Make a wide matrix with patient in columns 
prot_bl_tbl<-as.data.table(prot_bl_matrix)
# is there are duplicated this will return counts that are double 
check_dups<-data.table::dcast(prot_bl_tbl,  TESTNAME ~ PATNO,
                                value.var ='TESTVALUE')
subset(check_dups, select=-c(TESTNAME, PATNO))

hist(as.data.frame(subset(check_dups, -c(TESTNAME))))
prot_bl_wide<-data.table::dcast(prot_bl_tbl,  TESTNAME ~ PATNO,
                                value.var ='TESTVALUE', fun.aggregate = mean)


prot_bl_wide<-as.data.frame(prot_bl_wide)
rownames(prot_bl_wide)<-prot_bl_wide$TESTNAME
prot_bl_wide$TESTNAME<-NULL
row.names(prot_bl_wide)


fwrite(prot_bl_wide, paste0(output_files, 'untargeted_prot_bl.csv'), row.names = TRUE)


