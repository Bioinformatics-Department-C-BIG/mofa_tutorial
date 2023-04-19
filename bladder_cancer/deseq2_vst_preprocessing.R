
library(edgeR)
library(limma)
library(Glimma)
BiocManager::install('Seurat')
library('Seurat')
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library(DESeq2)
library("vsn")

library("SummarizedExperiment")


library('Seurat')

if (!require("pacman")) install.packages("pacman")
#BiocManager::install("vsn")

#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source('bladder_cancer/preprocessing.R')
source('bladder_cancer/preprocessing.R')

## Output directory

output_1='bladder_cancer/plots/deseq/'
output_files<-'bladder_cancer/'
output_de=paste0(output_1, 'gene')


##### Load required data 
Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)

seqdata<- t(X1_t_cut) ; rownames(seqdata)<-colnames(X1_t_cut)
countdata<-seqdata
sample_info<-Y_raw


#### Remove low expression 
raw_counts<-countdata
## filterbyExpr takes cpm so remove from there 
#cpm_data<-cpm(countdata, log=FALSE)

idx <- edgeR::filterByExpr(raw_counts, group = sample_info$Group)
raw_counts <- raw_counts[idx, ]


##### Define
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


vsn::meanSdPlot(log10(assay(dds)))
vsn::meanSdPlot(log10(counts_only))
vsn::meanSdPlot(counts_only)


# Compute normalization factors and vst 

dds <- estimateSizeFactors(dds,)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function
vsd <- varianceStabilizingTransformation(dds)


vsd_mat <- assay(vsd)
colnames(vsd_mat)<-vsd$Sample
vsd_mat
meanSdPlot(vsd)


##### Checks
# Check the effect of vst before and after
# Check distributions of samples using boxplots
boxplot(log10(assay(dds)), xlab="", ylab="Log2 counts ",las=2)

# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log10(assay(dds))),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Check distributions of samples using boxplots
boxplot(vsd_mat, xlab="", ylab="vst(counts) ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(vsd_mat),col="blue")
title("Boxplots of logCPMs (after vst)")
dev.off()

## ASSESS BATCH
### Assess batch effect
plotPCA(vsd, "Sample") + labs(color='Sample') + ggtitle("Batch effect") 


##### Store the most variable genes only for MOFA 
# Select most variable genes
vsd_mat<-selectMostVariable(vsd_mat, ng_g/100)
dim(vsd_mat)
highly_variable_genes_mofa<-vsd_mat
  
write.csv(highly_variable_genes_mofa,'highly_variable_genes_mofa.csv')
write.csv(highly_variable_genes_mofa, paste0(output_files,'/highly_variable_genes_mofa.csv'))
write.table(t(vsd_mat), paste0(output_files,'/highly_variable_genes_mofa_t.txt'), row.names = FALSE, sep = '\t')
  

# Check that the distribution is approximately normal
hist(highly_variable_genes_mofa)



# TODO check seurat 
##### Seurat option

var_features<-Seurat::FindVariableFeatures(raw_counts, 
                                        selection.method='vst', nfeatures=2000)

var_features
top10 <- head(VariableFeatures(var_features), 10)


plot1 <- VariableFeaturePlot(var_features) + 
  theme(legend.position="top")

