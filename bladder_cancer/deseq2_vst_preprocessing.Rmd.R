
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

if (!require("pacman")) install.packages("pacman")
#BiocManager::install("vsn")

#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source('bladder_cancer/preprocessing.R')
source('bladder_cancer/preprocessing.R')


output_1='bladder_cancer/plots/deseq/'
output_files<-'bladder_cancer/'


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)

# Use this script just for deseq 1
prot=FALSE

#####

seqdata<- t(X1_t_cut) ; rownames(seqdata)<-colnames(X1_t_cut)
countdata<-seqdata
output_de=paste0(output_1, 'gene')





sample_info<-Y_raw


# remove low expression 
raw_counts<-countdata
NROW(raw_counts)

idx <- edgeR::filterByExpr(raw_counts[,1:ncol(raw_counts)], group = sample_info$Group)
raw_counts <- raw_counts[idx, ]

### batch effect and normalization 
# Create a separate matrix with counts only
counts_only <- raw_counts
# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)
dds <- DESeqDataSetFromMatrix(
  countData = counts_only,
  colData = sample_info,
  design = ~Sample, tidy = F  
  
)

NROW(dds)
log10(assay(dds))

vsn::meanSdPlot(log10(assay(dds)))
vsn::meanSdPlot(log10(counts_only))
vsn::meanSdPlot(counts_only)


# Compute normalization factors
dds <- estimateSizeFactors(dds,)
sizeFactors(dds)

# Variance stabilization transformation
vsd <- varianceStabilizingTransformation(dds)
vsd_mat <- assay(vsd)

vsd_mat

 #vsd <- vst(dds, nsub=nrow(dds))
# vsd_mat <- assay(vsd)

meanSdPlot(vsd)

# Check the effect of vst before and after
# Check distributions of samples using boxplots
boxplot(log(assay(dds)), xlab="", ylab="Log2 counts ",las=2)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log(assay(dds))),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Check distributions of samples using boxplots
boxplot(vsd_mat, xlab="", ylab="Log2 counts ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(vsd_mat),col="blue")
title("Boxplots of logCPMs (after vst)")

# 
top_perc=0.25
vsd_mat<-select_most_variable(vsd_mat, 0.25)
##### Select most variable genes


dev.off()



## ASSESS BATCH
### Assess batch effect

plotPCA(vsd, "Sample") + labs(color='Sample') + ggtitle("Batch effect") 

## save output
  highly_variable_genes_mofa<-vsd_mat
  write.csv(highly_variable_genes_mofa,'highly_variable_genes_mofa.csv')
  write.csv(highly_variable_genes_mofa, paste0('bladder_cancer/highly_variable_genes_mofa.csv'))
  write.table(t(vsd_mat), 'bladder_cancer/highly_variable_genes_mofa_t.txt', row.names = FALSE, sep = '\t')
  

# todo need to apply just vsn to this one
hist(highly_variable_proteins_mofa)
hist(highly_variable_genes_mofa)

