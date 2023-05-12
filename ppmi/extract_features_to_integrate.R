
# TODO: extract all visits at once, then separate them later on for your analysus
library(dplyr)
library(data.table)
library(stringr)
source(paste0(script_dir,'/setup_os.R'))
os_dir='os_dir'
#rnas<-read.csv2('ppmi/ppmi_data/rnaseq/featCounts_SL_1239.longRNA_20230217.csv', sep = ',')
#### too many files so we separate each time points 
##### Maybe do them one by one in a function



rnas<-read.csv2(paste0('ppmi/ppmi_data/rnaseq/', VISIT, '.csv'), sep = ',')

output_files='ppmi/output/'


VISIT='V08'
VISIT='BL'
Visits<-c('BL', 'V04', 'V06', 'V08')

getwd()
for (VISIT in Visits){
  rnas<-read.csv2(paste0('ppmi/ppmi_data/rnaseq/', VISIT, '.csv'), sep = ',')
  rownames(rnas)<-rnas$Geneid
  rnas$Geneid<-NULL
  rnas_BL<-select(rnas,contains(VISIT))
  
  ### Split the names to extract patient number
  names_split<- strsplit(names(rnas_BL),split='\\.')
  names_split_df<-do.call(rbind, names_split);names_split_df
  #PATNO<-sapply(names_split, '[[', 2)
  PATNO<-names_split_df[,2]
  length(unique(PATNO))
  
  dim(rnas_BL)
  colnames(rnas_BL)<-PATNO;head(colnames(rnas_BL))
  
  write.csv2(rnas_BL, paste0(output_files, 'rnas_', VISIT, '.csv'), row.names = TRUE)
  
}


# Also merge all RNA files  for all visits together 
#rna_all_visits<-sapply(visits, function(VISIT){
#  read.csv2(paste0('ppmi/ppmi_data/rnaseq/', VISIT, '.csv'),
#            sep = ',')}
#        )

rna_all_visits_list<-sapply(Visits, 
  ### here i just merged them                          
      function(VISIT){
    ### RNA visits all 
  rnas_file = paste0(output_files, 'rnas_', VISIT, '.csv')
    rnas<-as.matrix(fread(rnas_file, header=TRUE), 
              rownames=1)
    # PATNO_VISIT column index 
    colnames(rnas)<-paste0(colnames(rnas), '_', VISIT)
    return(rnas)
    
}
  
)
# bind all visits for all patients together with underscore PATNO_VISIT
rna_all_visits<-do.call(cbind,rna_all_visits_list )

write.csv2(rna_all_visits, paste0(output_files, 'rnas_all_visits.csv'), row.names = TRUE)


#### Read in the mirnas 
mirnas_rpmmm<-read.csv2('ppmi/ppmi_data/mirnas/PPMI_sncRNAcounts/mirna_quantification_matrix_rpmmm_norm.csv/std_quantification_rpmmm_norm_mirna.final_ids.csv', sep = '\t')
mirnas_rpmmm<-read.csv2('ppmi/ppmi_data/mirnas/PPMI_sncRNAcounts/mirna_quantification_matrix_raw.csv/std_quantification_raw_mirna.final_ids.csv', sep = '\t')

names<-colnames(mirnas_rpmmm)[-1]

#names_split<- strsplit(names,split='\\.')
#visit<-sapply(names_split, '[[', 4)
#bl_visits<-visit=='BL'

VISIT='V08'
VISIT='BL'
VISIT='V04'
VISIT='V06'
mirnas_BL<-select(mirnas_rpmmm,contains(VISIT))
mirnas_BL

names_split<- strsplit(names(mirnas_BL),split='\\.')

PATNO<-sapply(names_split, '[[', 3)
length(PATNO)
colnames(mirnas_BL)<-PATNO
rownames(mirnas_BL)<-mirnas_rpmmm$miRNA

length((mirnas_BL))
tail(colnames(mirnas_BL), 400)
rownames(mirnas_BL)


write.csv2(mirnas_BL, paste0(output_files, 'mirnas_', VISIT, '.csv'), row.names = TRUE)


paste0(output_files, 'mirnas_', VISIT, '.csv')

names<-colnames(mirnas_rpmmm)[-1]


######### ADD all visits
## TODO: remove the code above and make separate files from the bigger one 
VISIT='V08'
VISIT='BL'

mirnas_all_visits<-mirnas_rpmmm
mirnas_df<-mirnas_rpmmm[-1]
#cbind(mirnas_df[,1],mirnas_rpmmm[,2])
# split column names 


names_split<- strsplit(names(mirnas_df),split='\\.')
names_split_df<-do.call(rbind, names_split)

dim(names_split_df);length(colnames(mirnas_df))
PATNO<-names_split_df[,3]
VISIT<-names_split_df[,4]
PATNO_VISIT<-paste0(PATNO, '_', VISIT); length(PATNO_VISIT); length(unique(PATNO_VISIT))
length(unique(paste0(PATNO, '_', VISIT)))
colnames(mirnas_df)<-PATNO_VISIT
length(unique(colnames(mirnas_df))); length(PATNO_VISIT)
rownames(mirnas_df)<-mirnas_rpmmm$miRNA; head(rownames(mirnas_df)); head(mirnas_df)


### TODO: which are the duplicates? can we delete them? 
length((mirnas_df))
tail(colnames(mirnas_df))
rownames(mirnas_df)
dim(mirnas_df)
length(colnames(rna_all_visits)); 
length(unique(colnames(rna_all_visits)))
write.csv2(mirnas_df, paste0(output_files, 'mirnas_all_visits.csv'), row.names = TRUE)



library(tidyr)
##### PROTEOMICS - OLINK



TISSUE='CSF'
TISSUE='Plasma'

NORMALIZED=TRUE


setwd(os_dir)
Visits


output_files<-paste0(data_dir, 'ppmi/output/')
NORMALIZED=TRUE

Visits=c('BL', 'V04', 'V06', 'V08')

#### ALL VISITS # include alla visits again

  prot_files<-list.files(path=paste0(data_dir,'ppmi/ppmi_data/proteomics/targeted_olink/', TISSUE), pattern='*_NPX*',
                         full.names = TRUE)
  
  if (NORMALIZED){
    prot_files<-list.files(path=paste0(data_dir,'ppmi/ppmi_data/proteomics/targeted_olink/', TISSUE), pattern='*_NPX*',
                           full.names = TRUE)
    
    pv='NPX'
  }else{
    
    prot_files<-list.files(path=paste0('ppmi/ppmi_data/proteomics/targeted_olink/', TISSUE), pattern='*Counts*',
                           full.names = TRUE)
    pv='COUNT' # Value of the protein level 

    
  }
  
  #prot_files_1<-read.csv('ppmi/ppmi_data/proteomics/targeted_olink/Plasma/PPMI_Project_196_Plasma_Cardio_NPX.csv')
  #prot_files_2<-read.csv('ppmi/ppmi_data/proteomics/targeted_olink/Plasma/PPMI_Project_196_Plasma_Cardio_Counts.csv')
  
  #intersect(unique(prot_files_1$PATNO),unique(prot_files_2$PATNO))
  
  
  prot_files
  outname<-paste0(output_files, 'proteomics_', TISSUE, '_',NORMALIZED,  '.csv')
  outname2<-paste0(output_files, 'proteomics_', TISSUE, '_', NORMALIZED,  '_no_log.csv')
  
  
  # Merge the 4 datasets together 
  read_all<-lapply(prot_files, function(x) {
    # panel INF is read as a number...
    df<-read.csv(x)
    df$PANEL<-as.character(df$PANEL)
    return(df)
  }
                                )
  
  all_frames<-lapply(read_all, as.data.frame)
  all_frames2<-do.call(rbind, all_frames)

  ppmi_prot=all_frames2
  length(unique(all_frames2$PATNO))
  length(unique(all_frames2$ASSAY))
  ### remove PPMI- suffix from PATNO column 
  ppmi_prot$PATNO<-str_replace(ppmi_prot$PATNO,'PPMI-', '')
  print(paste0('Unique patients in data: ', length(unique(ppmi_prot$PATNO))))
  prot_bl<-ppmi_prot
  
  
  
  #ppmi_prot<-as.data.frame(ppmi_prot) %>% 
  #               mutate(across('PATNO', str_replace, 
  #                          'PPMI-', ''
  #                         ))
  ## Filter baseline
  ## Remove unecessary columns 
  prot_bl_matrix<-as.data.frame(subset(prot_bl,
                                       select=-c(LOD, update_stamp, OLINKID, UNIPROT, MISSINGFREQ,
                                                 PANEL, EVENT_ID, PLATEID, INDEX, PANEL_LOT_NR)))
  
  prot_bl_matrix<-as.data.frame(subset(prot_bl,
                                       select=-c( update_stamp, OLINKID, UNIPROT,
                                                  PANEL, PLATEID)))
  ## Extract QC pass only 
    prot_bl_matrix
  prot_bl_matrix_QC<-prot_bl_matrix[prot_bl_matrix$QC_WARNING=='PASS',]
  if (length(prot_bl_matrix_QC==0)){
    prot_bl_matrix_QC<-prot_bl_matrix
    
  }
  hist(log10(prot_bl_matrix[,pv]))
  
  
  
  ## Reshape: Make a wide matrix with patient in columns 
  prot_bl_tbl<-as.data.table(prot_bl_matrix)
  NROW(unique(prot_bl_tbl$PATNO))
  prot_bl_tbl$PATNO_EVENT_ID<-paste0(prot_bl_tbl$PATNO,'_',prot_bl_tbl$EVENT_ID)

  prot_bl_wide<-data.table::dcast(prot_bl_tbl,  ASSAY ~ PATNO_EVENT_ID,
                                  value.var =pv, fun.aggregate = mean)
  
  length(colnames(prot_bl_wide)); length(unique(colnames(prot_bl_wide)))
  
  prot_bl_wide<-as.data.frame(prot_bl_wide)
  rownames(prot_bl_wide)<-prot_bl_wide$ASSAY
  prot_bl_wide$ASSAY<-NULL
  head(row.names(prot_bl_wide))
  
  
  ## Remove the log normalization
  prot_bl_wide_unlog<-as.data.frame(sapply(prot_bl_wide, function(x){2**(as.numeric(x))}),
                                    row.names =row.names(prot_bl_wide) )
  
  
  ## Write output both log and not logged
  fwrite(prot_bl_wide, paste0(outname), row.names = TRUE)
  fwrite(prot_bl_wide_unlog, paste0(outname2), row.names = TRUE)
  dim(prot_bl_wide)
  dim(prot_bl_wide)
  
  df<-prot_bl_wide_unlog[prot_bl_wide_unlog<10]
  #df<-prot_bl_wide[prot_bl_wide<200]
  
  hist(as.numeric(as.matrix(df)))

dev.off()


combined[match(colnames(prot_bl_wide), combined$PATNO_EVENT_ID),]$COHORT_DEFINITION








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
#check_dups<-data.table::dcast(prot_bl_tbl,  TESTNAME ~ PATNO,
  #                              value.var ='TESTVALUE')
#subset(check_dups, select=-c(TESTNAME, PATNO))

#hist(as.data.frame(subset(check_dups, -c(TESTNAME))))
prot_bl_wide<-data.table::dcast(prot_bl_tbl,  TESTNAME ~ PATNO,
                                value.var ='TESTVALUE', fun.aggregate = mean)


prot_bl_wide<-as.data.frame(prot_bl_wide)
rownames(prot_bl_wide)<-prot_bl_wide$TESTNAME
prot_bl_wide$TESTNAME<-NULL
row.names(prot_bl_wide)


fwrite(prot_bl_wide, paste0(output_files, 'untargeted_prot_bl.csv'), row.names = TRUE)


