
#BiocManager::install('vsn')
#install.packages('gplots')


library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library(DESeq2)

library("SummarizedExperiment")

library(data.table)
library(dplyr)
library(stringr)
## Output directory

output_1='ppmi/output/'
output_files_orig<-'ppmi/output/'
output_files<-'ppmi/output/'
outdir_orig<-'ppmi/plots/'


sanity_check_metadata_align<-function(counts, sample_info){
  
  cnames<-colnames(counts)
  snames_1<-str_split_fixed(cnames, '_',2)[,1]
  snames_2<-sample_info$Sample
  return( all(snames_1==snames_2))
  
}




output_de=paste0(output_1, 'gene')
#source('bladder_cancer/preprocessing.R')
# TODO: move the pre-processing script to utils


##### Load required data 
# TODO: input all the visits 

MIN_COUNT_G=100
MIN_COUNT_M=10
TOP_GN=0.10
TOP_MN=0.50

VISIT='BL'
VISIT='V04'
VISIT='V08'

g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(TOP_MN, '_', MIN_COUNT_M, '_') 




metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)

#### Remove low expression 
process_mirnas<-FALSE
if (process_mirnas){
  VISIT='BL'
  mirnas_file<-paste0(output_files, 'mirnas_',VISIT,  '.csv')
  mirnas_BL<-as.data.frame( as.matrix(fread(mirnas_file, header=TRUE), rownames=1))
  
  VISIT='V08'
  mirnas_file<-paste0(output_files, 'mirnas_',VISIT,  '.csv')
  mirnas_V08<-as.data.frame(as.matrix(fread(mirnas_file, header=TRUE), rownames=1))
  
  raw_counts<-as.data.frame(mirnas_BL)
  
  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  param_str_m<-paste0('mirnas_', m_params)
  highly_variable_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')
  
  
}else{
  
  rnas_file<-paste0(output_files, 'rnas_all_visits.csv')
  rnas_all_visits_r<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)

  
  
  # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
  #raw_counts<-raw_counts %>% select(common_samples)
  
  min.count=MIN_COUNT_G
  most_var=TOP_GN
  param_str_g<-paste0('rnas_', g_params )
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  outdir_s<-paste0(outdir_orig, '/single/', param_str_g, 'all_visits')

  
  
}

#### TODO: 
# BIND ALL VISITS AND TIME POINTS TOO!! 
# ONLY BIND WHATS common




##### Define

### TODO: Question: Should I input everything into the matrix to normalize? 
### And then filter 

### batch effect and normalization 
# Create a separate matrix with counts only



Sample_info_all<-str_split_fixed( colnames(rnas_all_visits_r), '_', 2);Sample=Sample_info_all[,1];Visit=Sample_info_all[,2]





### VST? The variance stabilizing and rlog transformations are provided for applications other than differential testing,
### for example clustering of samples or other machine learning applications. For differential testing we recommend the DESeq function applied to raw counts as outlined above.



#### NUMBER 2 COMPARE WITH CONTROLS 
##### separate cohorts and create filters 

# select patients 


df<-rnas_all_visits_r
dim(Sample_info_all)


### first select samples in all visits 
### Experiment 1: 





#### First filter the metadata and then order them!! 

### filtery by COHORTS 
combined_sel<-combined[combined$COHORT %in% c(1,2),]
dim(df); dim(combined_sel)
common_in_meta<-unique(intersect( Sample, combined_sel$PATNO ))




### Second filter the dataframe
# which of the common are in the counts matrix
length(Sample %in% common_in_meta) 



df_filt<-df[,(Sample %in% common_in_meta) ] # select only columns with common patients 
sample_info<-str_split_fixed( colnames(df_filt), '_', 2);Sample=sample_info[,1];Visit=sample_info[,2]
dim(sample_info)
dim(df_filt)


#which(combined_sel_uniq$PATNO_EVENT_ID=='10874_BL')
#which(combined_sel$PATNO=='10874')
#View(combined[which(combined$PATNO=='10874'),'COHORT'])


### And then order the metadata based on the df 
combined_sel$PATNO_EVENT_ID = paste0(combined_sel$PATNO, '_', combined_sel$EVENT_ID) 
dim(df_filt)
dim(combined_sel)



### now select metadata from samplesheet and order it as in dataframe
#which(is.na(combined_ordered$PATNO_EVENT_ID))
# if there are two rows take the second
# TODO: better to take the one with less NAs 
combined_sel_uniq<-combined_sel[!duplicated(combined_sel$PATNO_EVENT_ID,fromLast=TRUE),]
which(is.na(combined_sel_uniq$COHORT))
which(is.na(combined_sel_uniq$PATNO_EVENT_ID))

combined_ordered<-combined_sel_uniq[match(colnames(df_filt) ,  combined_sel_uniq$PATNO_EVENT_ID ),]
colnames(df_filt[,which(is.na(combined_ordered$COHORT))])

# sanity check
which(is.na(combined_sel_uniq$PATNO_EVENT_ID))
which(is.na(combined_ordered$PATNO_EVENT_ID))


dim(df_filt)
dim(combined_ordered)



sample_info<-str_split_fixed( colnames(df_filt), '_', 2);Sample=sample_info[,1];Visit=sample_info[,2]
sample_info$COHORT=as.factor(combined_ordered$COHORT)
sample_info$Visit=as.factor(Visit)


####### AFTER ALIGNMENT WITH METADATA, filter by expression , input group
#raw_counts<-as.data.frame(tps_merged)
raw_counts=df_filt
idx <- edgeR::filterByExpr(raw_counts, group=c(sample_info$COHORT,sample_info$VISIT ))
length(which(idx))
raw_counts_filt <- as.data.frame(as.matrix(raw_counts[idx, ]))
dim(raw_counts)

df_filt=raw_counts_filt



## SANITY CHECK
dim(df_filt)[2]==length(Sample)


#### Now get new sample info again




### Now define what will go inside deseq

non_na_ids<-which(!is.na(combined_ordered$PATNO_EVENT_ID))
combined_ordered<-combined_ordered[non_na_ids,]
df_filt_2<-df_filt[, non_na_ids]


counts_only<-df_filt_2
sample_info_split<-str_split_fixed( colnames(df_filt_2), '_', 2)
Sample=sample_info_split[,1]
Visit=sample_info_split[,2]

#### common pipeline again
sample_info<-DataFrame(Sample=Sample, Visit=Visit)
sample_info$COHORT=as.factor(combined_ordered$COHORT)


sample_info$AGE_AT_VISIT = as.factor(combined_ordered$AGE_AT_VISIT)

#View(cbind(combined_ordered_2$PATNO,colnames(df_filt)) )




#### 
## nextflowrun
counts_only_nf<-counts_only
counts_only_nf$gene_id<-rownames(counts_only_nf)
rownames(counts_only_nf)<-NULL

profile_rna_file<-paste0(output_files_orig, '/single/rnas.csv')
#write.csv(counts_only_nf,profile_rna_file)


#### 

samplesheet_file<-paste0(output_files_orig, '/single/samplesheet.csv')

write.csv(sample_info,samplesheet_file )
contrasts_file<-paste0(output_files_orig, '/single/contrasts.csv')
contrasts<-'COHORT,Visit'
write.csv(contrasts, contrasts_file )




#### Make appropriate visit filters 
VISIT='V06'
Visits=c('BL', 'V04', 'V06', 'V08')
paired<-FALSE


VISIT='V08'

for (VISIT in Visits){
  

  outdir_s<-paste0(outdir_orig, '/single/', param_str_g, 'all_visits_', VISIT,'_', paired)
  
  dim(counts_only)

  to_sel<-which(sample_info$Visit==VISIT)
  sample_info_final<-sample_info[to_sel,]
  counts_only_final=counts_only[to_sel]
  dim(sample_info)
  
  sample_info_final$COHORT=as.factor(sample_info_final$COHORT)
  sanity_check_metadata_align(counts_only_final, sample_info_final)
  
  
  outdir_s
  ###
  # TODO: assign the groups 
  # counts_only
  #dds <- DESeqDataSetFromMatrix(
  #  countData = round(counts_only_final),
  #  colData = sample_info_final,
  #  design = ~COHORT+COHORT:Visit, tidy = F
  #)
  

  ## in specific visits only compare cohorts 
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_only_final),
    colData = sample_info_final,
    design = ~COHORT, tidy = F
  )
  
  
  
  
  
  #### Run DE 
  deseq2Data <- DESeq(dds)
  dir.create(outdir_s)
  saveRDS(deseq2Data, paste0(outdir_s, '/deseq_results.RDS'))
}








#### Make appropriate visit filters 
VISIT='V06'
Visits=c('BL', 'V04', 'V06', 'V08')
paired<-TRUE

View(rbind(c(sample_info_final$Sample), colnames(counts_only_final)))





#### keep all visits together 

intersect()






