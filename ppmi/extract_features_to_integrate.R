

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




##### PROTEOMICS
prot_files<-list.files(path='ppmi/ppmi_data/proteomics/targeted_olink/plasma/', pattern='*_NPX*',
                       full.names = TRUE)




read_all<-lapply(prot_files, read.csv)

all_frames<-lapply(read_all, as.data.frame)
all_frames2<-bind_rows(all_frames)
ppmi_prot=all_frames2

### remove ppmi from PATNO column 
ppmi_prot<-ppmi_prot %>% mutate(across('PATNO', str_replace, 
                            'PPMI-', ''
                            ))

prot_bl<-ppmi_prot[ppmi_prot$EVENT_ID=='BL',]



library(tidyr)
prot_bl_matrix<-as.data.frame(subset(prot_bl,
                    select=-c(LOD, update_stamp, OLINKID, UNIPROT, MISSINGFREQ,
                              PANEL, EVENT_ID, PLATEID, INDEX, PANEL_LOT_NR)))

prot_bl_matrix<-prot_bl_matrix[prot_bl_matrix$QC_WARNING=='PASS',]

prot_bl_2<-spread(prot_bl_matrix, ASSAY, NPX)
length(unique(prot_bl_matrix$PATNO))
length(unique(prot_bl_matrix$PANEL_LOT_NR))


prot_bl_2<-as.data.frame(subset(prot_bl_2,
                        select=-c(EVENT_ID, INDEX, PLATEID, 
                                  QC_WARNING)))

hist(prot_bl_matrix$NPX)


## make a wide matrix with patient in columns 
prot_bl_tbl<-as.data.table(prot_bl_matrix)
prot_bl_wide<-data.table::dcast(prot_bl_tbl,  ASSAY ~ PATNO,
                         value.var ='NPX', fun.aggregate = mean)


rownames(prot_bl_wide)<-prot_bl_wide$ASSAY
prot_bl_wide$ASSAY<-NULL
row.names(prot_bl_wide)

fwrite(prot_bl_wide, paste0(output_files, 'proteomics_bl.csv'), row.names = TRUE)


#prot_bl_2_t<-t(prot_bl_2)p
# MAKE THE FIRST ROW colname
#colnames(prot_bl_2_t)<-prot_bl_2_t[1,]
#prot_bl_2_t<-prot_bl_2_t[-1,]
  ## TODO: drop worst duplicates with most missing values  







