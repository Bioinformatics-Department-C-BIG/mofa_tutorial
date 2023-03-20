

#### This script performs filter by expression, size factor estimation and vsn on the whole matrix of 
#### RNAs or miRNAS 
#### parameters: MIN_COUNT_M, MIN_COUNT_G define the parameters for filtering out genes with low counts 
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library(DESeq2)
library("vsn")

library("SummarizedExperiment")

library(data.table)
library(dplyr)
## Output directory

output_1='ppmi/output/'
output_files_orig<-'ppmi/output/'

output_de=paste0(output_1, 'gene')
source(paste0(script_dir, '/../bladder_cancer/preprocessing.R'))


metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)
# TODO: move the pre-processing script to utils


##### Load required data 
# TODO: input all the visits 



MIN_COUNT_G=100
MIN_COUNT_M=10
VISIT='BL'
VISIT='V04'



sel_coh=c(1)
sel_coh_s<-paste(sel_coh,sep='_',collapse='-')
sel_coh_s
VISIT_S=paste(VISIT,sep='_',collapse='-')

g_params<-paste0(VISIT_S, '_', TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0( VISIT_S, '_', TOP_MN, '_', MIN_COUNT_M, '_') 


#### Remove low expression 
process_mirnas<-TRUE
if (process_mirnas){
   mirnas_file<-paste0(output_files, 'mirnas_all_visits.csv')
   mirnas_BL<-as.matrix(fread(mirnas_file, header=TRUE), rownames=1)
   
   raw_counts<-mirnas_BL

  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  param_str_m<-paste0('mirnas_', m_params ,sel_coh_s, '_')
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, 'mirnas_', param_str_m, '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')
  
  
}else{
  rnas_file<-paste0(output_files, 'rnas_all_visits.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
 
  raw_counts<-rnas_BL
    # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
# raw_counts<-raw_counts %>% select(common_samples)

  min.count=MIN_COUNT_G
  most_var=TOP_GN
  param_str_g<-paste0('rnas_', g_params, sel_coh_s, '_'  )
  vsn_out_file<-highly_variable_outfile<-paste0(output_files, 'rnas_', param_str_g,  '_vsn.csv')
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  
  highly_variable_outfile
  
}

rownames(raw_counts)

##### 1.  First create the summarized experiment object  
### find common samples in mirnas file + metadata
## subset and order by common samples
## And create SE object with metadata
duplicate_samples<-colnames(mirnas_BL)[which(duplicated(colnames(mirnas_BL),fromLast = TRUE))]

raw_counts_all=raw_counts
class(raw_counts_all) <- "numeric"
## They seem to have taken averages for replicas so need to fix 
raw_counts_all<-round(raw_counts_all)


## Question: why are there duplicate samples - seems to be controls! 
## first filter what is in metadata and mirnas ?


### TODO: move to a utils / preprocessing file because it is used also for proteoomics
library(SummarizedExperiment)
getSummarizedExperimentFromAllVisits<-function(raw_counts_all, combined){
  #
  raw_counts_all<-raw_counts_all[,!duplicated(colnames(raw_counts_all), fromLast=TRUE)]
  combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)
  
  ### some samples do not exist in metadata so filter them out 
  ## 
  common_samples<-intersect(colnames(raw_counts_all),combined$PATNO_EVENT_ID)
  unique_s<-colnames(raw_counts_all)[!(colnames(raw_counts_all) %in% common_samples)]
  metadata_filt<-combined[match(common_samples, combined$PATNO_EVENT_ID),]
  raw_counts_filt<-raw_counts_all[,match(common_samples, colnames(raw_counts_all))]
  dim(metadata_filt)[1] ==dim(raw_counts_filt)[2]
  
  
  #subset sample names
  raw_counts<-raw_counts_filt
  
  se=SummarizedExperiment(raw_counts_filt, colData = metadata_filt)
  return(se)
}


se<-getSummarizedExperimentFromAllVisits(raw_counts_all, combined)
# remove duplicates 


##### Up till here it is generic, no filters yet. 




##### 2.  Now start filtering to normalize as appropriate 
## Option 1: normalize cohort and EVENT separately!! 



se_filt<-se[,(se$EVENT_ID %in% VISIT & se$COHORT %in% sel_coh )]

Sample<-colnames(se_filt)
sample_info<-DataFrame(Sample=Sample)

raw_counts=assays(se_filt)[[1]]

## filterbyExpr takes cpm so remove from there 
idx <- edgeR::filterByExpr(raw_counts,min.count=min.count)

length(which(idx))
raw_counts <- as.matrix(raw_counts[idx, ])
dim(raw_counts)
se_filt=se_filt[idx]
se_filt


##### Define

### TODO: Question: Should I input everything into the matrix to normalize? 
### And then filter 

### batch effect and normalization 
# Create a separate matrix with counts only
# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)


### DEFINE THE DESEQ OBJECT with the groups appropriately 
se_filt$EVENT_ID=as.factor(se_filt$EVENT_ID)
se_filt$COHORT=as.factor(se_filt$COHORT)
se_filt$PATNO=as.factor(se_filt$PATNO)

if (length(sel_coh)>1){
  
      if (length(VISIT)>1){
        ddsSE <- DESeqDataSet(se_filt, 
                              design = ~COHORT + EVENT_ID)
        vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
        
      }else{
      ddsSE <- DESeqDataSet(se_filt, 
                            design = ~COHORT)
      vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
    
  }
  }else{
    ddsSE <- DESeqDataSet(se_filt, 
                          design = ~PATNO)
    vsd <- varianceStabilizingTransformation(ddsSE)
    
  }




# Compute normalization factors and vst 
# or use blind=false 

#ddsSE <- estimateSizeFactors(ddsSE)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function


vsd_mat <- assay(vsd)

meanSdPlot(vsd_mat)


##### Checks
# Check the effect of vst before and after
par(mfrow=c(1,3))
# Check distributions of samples using boxplots
boxplot(log2(assay(ddsSE)[,1:30]), xlab="", ylab="Log2 counts ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log10(assay(dds))),col="blue")
title("Boxplots of logCPMs (unnormalised)")
boxplot(log10(raw_counts)[,1:30], xlab="", ylab="Log10 counts ",las=2)

# Check distributions of samples using boxplots
boxplot(vsd_mat[,1:30], xlab="", ylab="vst(counts) ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(vsd_mat),col="blue")
title("Boxplots of logCPMs (after vst)")

## ASSESS BATCH
### Assess batch effect
#plotPCA(vsd, "sample") + labs(color='sample') + ggtitle("Batch effect") 
#dev.off()

#### This saves the whole file without filtering for highly variable 
write.csv(vsd_mat,vsn_out_file)
##### Store the most variable genes only for MOFA 
# Select most variable genes
### run on its own for all visits? 
# Check that the distribution is approximately normal
dev.off()

###TODO: Move this to the mofa file 
highly_variable_genes_mofa<-selectMostVariable(vsd_mat, most_var)
write.csv(highly_variable_genes_mofa, highly_variable_outfile, col.names = TRUE)
dim(highly_variable_genes_mofa)
rownames(highly_variable_genes_mofa)


### SANITY CHECK: Just plot one gene before and after preprocessing to ensure the mapping looks correct 
df=highly_variable_genes_mofa
par(mfrow=c(2,1))
idx=30
plot(df[idx,1:150]) 
gname<-rownames(df)[idx]
title(gname)
df=raw_counts
plot(df[gname, 1:150])
title(gname)

