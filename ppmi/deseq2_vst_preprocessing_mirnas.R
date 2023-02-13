
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


## Output directory

output_1='ppmi/output/'
output_files<-'ppmi/output/'
output_de=paste0(output_1, 'gene')
source('bladder_cancer/preprocessing.R')
# TODO: move the preprocessing script to utils


##### Load required data 
# TODO: input all the visits 



#### Remove low expression 


raw_counts<-as.data.frame(mirnas_BL)
Sample<-colnames(mirnas_BL)
sample_info<-DataFrame(Sample=Sample)

raw_counts<-mutate_all(raw_counts, as.numeric)


hist(log2(as.matrix(raw_counts)))

## filterbyExpr takes cpm so remove from there 
#cpm_data<-cpm(countdata, log=FALSE)

idx <- edgeR::filterByExpr(raw_counts,)
raw_counts <- as.matrix(raw_counts[idx, ])


##### Define

### TODO: Question: Should I input everything into the matrix to normalize? 
### And then filter 

### batch effect and normalization 
# Create a separate matrix with counts only
counts_only <- raw_counts
# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_only),
  colData = sample_info,
  design = ~Sample, tidy = F
)




# Compute normalization factors and vst 

dds <- estimateSizeFactors(dds,)
sizeFactors(dds)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function
vsd <- varianceStabilizingTransformation(dds)


vsd_mat <- assay(vsd)
colnames(vsd_mat)<-vsd$Sample
vsd_mat
meanSdPlot(vsd_mat)


##### Checks
# Check the effect of vst before and after
# Check distributions of samples using boxplots
boxplot(log10(assay(dds)), xlab="", ylab="Log2 counts ",las=2)
dev.off()
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log10(assay(dds))),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Check distributions of samples using boxplots
boxplot(vsd_mat[,1:30], xlab="", ylab="vst(counts) ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(vsd_mat),col="blue")
title("Boxplots of logCPMs (after vst)")


## ASSESS BATCH
### Assess batch effect
plotPCA(vsd, "Sample") + labs(color='Sample') + ggtitle("Batch effect") 


##### Store the most variable genes only for MOFA 
# Select most variable genes
highly_variable_genes_mofa<-selectMostVariable(vsd_mat, 0.9)
dim(highly_variable_genes_mofa)


write.csv(highly_variable_genes_mofa,'highly_variable_genes_mofa.csv')
write.csv(highly_variable_genes_mofa, paste0(output_files,'/highly_variable_genes_mofa.csv'), col.names = TRUE)
write.table(t(vsd_mat), paste0(output_files,'/highly_variable_genes_mofa_t.txt'), row.names = FALSE, sep = '\t')

dim(highly_variable_genes_mofa)

# Check that the distribution is approximately normal
hist(highly_variable_genes_mofa)


