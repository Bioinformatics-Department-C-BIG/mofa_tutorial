




library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)

if (!require("pacman")) install.packages("pacman")
#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source('bladder_cancer/preprocessing.R')


output_1='bladder_cancer/plots/deseq/'
output_files<-'bladder_cancer/'


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)

prot=FALSE


if (prot){
  lowest_thresh<-10
  params.remove_low_threshold<-15
  params.ng<-2
}else{
  lowest_thresh<-10
  params.remove_low_threshold<-10
  params.ng<-5
  
  
}



if (prot){
  
  seqdata <- read.delim(paste0(dir,'Proteomics_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X2_t_cut) ; rownames(seqdata)<-colnames(X2_t_cut)
  no_name<-which(is.na(seqdata[1]))
  
  
  #### Format
  if (length(no_name)>0){
    seqdata<-seqdata[,-no_name]}
  countdata<-seqdata
  
  output_de=paste0(output_1, 'prot')
  
  
}else{
  seqdata <- read.delim(paste0(dir,'RNAseq_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X1_t_cut) ; rownames(seqdata)<-colnames(X1_t_cut)
  no_name<-which(is.na(seqdata[1]))
  #### Format
  if (length(no_name)>0){
    seqdata<-seqdata[,-no_name]}
  countdata<-seqdata
  
  
  
  output_de=paste0(output_1, 'genes')
  
}


dim(seqdata)


seqdata
sample_info<-Y_raw


# remove low expression 
raw_counts<-countdata
idx <- edgeR::filterByExpr(raw_counts[,1:ncol(raw_counts)], group = sample_info$Group)
raw_counts <- raw_counts[idx, ]
# remove low variance
variances <- apply(raw_counts, 1, var)
topx<-names(variances[order(variances, decreasing = TRUE)])[1:(length(variances)/4)]
raw_counts <- raw_counts[topx, ]
NROW(raw_counts)

hist(variances)
dev.off()
### batch effect and normalization 
# Create a separate matrix with counts only
counts_only <- raw_counts
# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)
dds <- DESeqDataSetFromMatrix(
  countData = round( counts_only, digits=0),
  colData = sample_info,
  design = ~Sample, tidy = F  
  
)
library("vsn")
meanSdPlot(assay(ntd))

NROW(dds)
# Compute normalization factors
dds <- estimateSizeFactors(dds,)
sizeFactors(dds)


vsd <- vst(dds)
vsd_mat <- assay(vsd)
 

vsd <- vst(dds, nsub=nrow(dds))
vsd_mat <- assay(vsd)


# Apply VST normalization
#meanSdPlot(assay(vsd))

# Check distributions of samples using boxplots
boxplot(log(assay(dds)), xlab="", ylab="Log2 counts ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log(assay(dds))),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Check distributions of samples using boxplots
boxplot(vsd_mat, xlab="", ylab="Log2 counts ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(vsd_mat),col="blue")
title("Boxplots of logCPMs (unnormalised)")





## ASSES BATCH
### Assess batch effect

plotPCA(vsd, "Sample") + labs(color='Sample') + ggtitle("Batch effect") 





if (prot){
  highly_variable_proteins_mofa<-vsd_mat
  write.csv(highly_variable_proteins_mofa,'highly_variable_proteins_mofa.csv')
}else{
  highly_variable_genes_mofa<-vsd_mat
  write.csv(highly_variable_genes_mofa,'highly_variable_genes_mofa.csv')
}

hist(highly_variable_proteins_mofa)
hist(highly_variable_genes_mofa)

NROW(highly_variable_proteins_mofa)

