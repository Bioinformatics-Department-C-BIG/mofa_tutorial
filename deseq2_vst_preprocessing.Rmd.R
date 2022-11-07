
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
idx <- edgeR::filterByExpr(raw_counts[,1:ncol(raw_counts)], group = sample_info$Group)
raw_counts <- raw_counts[idx, ]

### batch effect and normalization 
# Create a separate matrix with counts only
counts_only <- raw_counts
# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_only, digits=0),
  colData = sample_info,
  design = ~Sample, tidy = F  
  
)

NROW(dds)
# Compute normalization factors
dds <- estimateSizeFactors(dds,)
sizeFactors(dds)

# Variance stabilization transformation
vsd <- vst(dds)
vsd_mat <- assay(vsd)

 

vsd <- vst(dds, nsub=nrow(dds))
vsd_mat <- assay(vsd)

# Check the effect of vsd before and after
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


##### Select most variable genes
variances <- apply(vsd_mat, 1, var)
topx<-names(variances[order(variances, decreasing = TRUE)])[1:(length(variances)/4)]
vsd_mat <- vsd_mat[topx, ]
NROW(vsd_mat)

hist(variances)
dev.off()



## ASSESS BATCH
### Assess batch effect

plotPCA(vsd, "Sample") + labs(color='Sample') + ggtitle("Batch effect") 

## save output
  highly_variable_genes_mofa<-vsd_mat
  write.csv(highly_variable_genes_mofa,'highly_variable_genes_mofa.csv')


# todo need to apply just vsn to this one
hist(highly_variable_proteins_mofa)
hist(highly_variable_genes_mofa)

