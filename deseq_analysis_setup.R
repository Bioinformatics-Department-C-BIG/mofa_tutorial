
script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))

print(script_dir)
library('R.filesets')
library(DESeq2)
library("SummarizedExperiment")
library(data.table)
library(dplyr)


### TODO: Add volcano plot for each time point -DONE
### TODO: add heatmap for all tps tpogether -DONE
#source('ppmi/de')

#load
library("factoextra")
library("FactoMineR")
library('pheatmap')
library('ggplot2')

#### Run DE 


#ddsSE <- DESeqDataSet(se_filt, 
#                      design = ~PATNO)

# TODO: assign the groups 
#dds <- DESeqDataSetFromMatrix(
# countData = assay(se_filt),
#  colData = colData(se_filt),
#  design = ~COHORT, tidy = F
#)

#source(paste0(script_dir, '/config.R'))

### LOAD runs


source(paste0(script_dir, '/../bladder_cancer/preprocessing.R'))
source(paste0(script_dir, '/utils.R'))

VISIT='V08'


process_mirnas<-FALSE


source(paste0(script_dir, '/config.R'))



print(deseq_file)

datalist=loadRDS(deseq_file)
ddsSE=datalist[[1]]
vsd=datalist[[2]]
se_filt=datalist[[3]]
deseq2Results=datalist[[4]]

table(se_filt$COHORT_DEFINITION)

# todo join strings
# TODO: Report the number of samples too! 

des<-gsub(' ', '', paste0(as.character(design(ddsSE))[-1]))

#if  (process_mirnas){
#
#  outdir_s<-paste0(outdir_orig, '/single/', param_str_m_f, des)
#  
#}else{
#  outdir_s<-paste0(outdir_orig, '/single/', param_str_g_f, des)
#  
#}
#outdir_s
# RUN DIFFERENTIAL EXPRESSION ANALYSIS 


### TODO: save the file so we don't have to fit the model each time!! 
dds<-ddsSE



#rm(deseq2Data)



dir.create(outdir_s)
#deseq2Data<-loadRDS(paste0(outdir_s, '/deseq_results.RDS'))
#### First obtain the single omics significant RNAs 



write.csv(deseq2Results, paste0(outdir_s, '/results.csv'))

deseq2ResDF <- as.data.frame(deseq2Results)
deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
deseq2ResDF$abslog2pval<-abs(deseq2ResDF$log2pval)

write.csv(deseq2ResDF, paste0(outdir_s, '/results_df.csv'))


### Up to here output can be used for other reasons
##

