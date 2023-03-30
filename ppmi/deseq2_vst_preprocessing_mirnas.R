
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
source('bladder_cancer/preprocessing.R')
# TODO: move the pre-processing script to utils


##### Load required data 
# TODO: input all the visits 

MIN_COUNT_G=100
MIN_COUNT_M=10
TOP_GN=0.10
TOP_MN=0.50

VISIT='BL'
VISIT='V08'
sel_coh='1'




g_params<-paste0(TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0( TOP_MN, '_', MIN_COUNT_M, '_') 


#### Remove low expression 
process_mirnas<-FALSE
if (process_mirnas){
   mirnas_file<-paste0(output_files, 'mirnas_',VISIT,  '.csv')
   mirnas_BL<-as.matrix(fread(mirnas_file, header=TRUE), rownames=1)
  
    raw_counts<-as.data.frame(mirnas_BL)

  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  param_str_m<-paste0('mirnas_', m_params)
  highly_variable_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')
  
  
}else{
  rnas_file<-paste0(output_files, 'rnas_all_visits.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
 
  raw_counts<-as.data.frame(rnas_BL)
    # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
# raw_counts<-raw_counts %>% select(common_samples)

  min.count=MIN_COUNT_G
  most_var=TOP_GN
  param_str_g<-paste0('rnas_', g_params )
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  
  
  
}


Sample<-colnames(raw_counts)
sample_info<-DataFrame(Sample=Sample)

raw_counts<-mutate_all(raw_counts, as.numeric)

## filterbyExpr takes cpm so remove from there 
#cpm_data<-cpm(countdata, log=FALSE)
idx <- edgeR::filterByExpr(raw_counts,min.count=min.count)

length(which(idx))
raw_counts <- as.matrix(raw_counts[idx, ])
dim(raw_counts)


raw_counts
##### Define

### TODO: Question: Should I input everything into the matrix to normalize? 
### And then filter 

### batch effect and normalization 
# Create a separate matrix with counts only
counts_only <- raw_counts

# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)

# TODO: assign the groups 
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_only),
  colData = sample_info,
  design = ~Sample, tidy = F
)




# Compute normalization factors and vst 

dds <- estimateSizeFactors(dds,)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function
vsd <- varianceStabilizingTransformation(dds)


vsd_mat <- assay(vsd)
colnames(vsd_mat)<-vsd$Sample

meanSdPlot(vsd)


##### Checks
# Check the effect of vst before and after
par(mfrow=c(1,3))
# Check distributions of samples using boxplots
boxplot(log2(assay(dds)[,1:30]), xlab="", ylab="Log2 counts ",las=2)
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

##### Store the most variable genes only for MOFA 
# Select most variable genes
highly_variable_genes_mofa<-selectMostVariable(vsd_mat, most_var)
dim(highly_variable_genes_mofa)


boxplot(highly_variable_genes_mofa[,1:30], xlab="", ylab="vst(counts) ",las=2)


write.csv(highly_variable_genes_mofa, highly_variable_outfile, col.names = TRUE)


dim(highly_variable_genes_mofa)

# Check that the distribution is approximately normal
dev.off()

hist(highly_variable_genes_mofa)




### SANITY CHECK: Just plot one gene before and after to ensure the mapping looks correct 
df=highly_variable_genes_mofa
idx=30
plot(df[idx,1:20]) 
gname<-rownames(df)[idx]
title(gname)
df=raw_counts
plot(df[gname, 1:20])
title(gname)

