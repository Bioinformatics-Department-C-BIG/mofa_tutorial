
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library("vsn")
library("DEP")
#BiocManager::install('DEP')

if (!require("pacman")) install.packages("pacman")
#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source('bladder_cancer/preprocessing.R')


output_1='bladder_cancer/plots/deseq/'
output_files<-'bladder_cancer/'


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)

prot=TRUE
# TODO: download LIBRARY DEP
# TODO: make a rmarkdown!

  seqdata <- read.delim(paste0(dir,'Proteomics_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X2_t_cut) ; rownames(seqdata)<-colnames(X2_t_cut)
  no_name<-which(is.na(seqdata[1]))
  
  
 
  countdata<-seqdata
  
  output_de=paste0(output_1, 'prot')
  #


dim(seqdata)


seqdata
sample_info<-Y_raw


# remove low expression 
raw_counts<-countdata
idx <- edgeR::filterByExpr(raw_counts[,1:ncol(raw_counts)], group = sample_info$Group)
raw_counts <- raw_counts[idx, ]



data<-raw_counts
data<-as.data.frame(data)
colnames(data)<-Y_raw$Sample

data$name<-c(rownames(data))
data$ID<-data$name

data$name
data$ID
data$ID# Make SummarizedExperiment
columns=seq(1:16)
exp_design = data.frame(label=Y_raw$Sample,condition=Y_raw$Sample, replicate=rep(1, 16))
se <- make_se(data, columns,exp_design)

# Filter and normalize


normalized_data<-normalize_vsn(se)
 meanSdPlot(normalized_data)
 

boxplot(normalized_data)

vsn_mat<-assays(normalized_data)[[1]]


##### Select most variable genes
variances <- apply(vsn_mat, 1, var)
topx<-names(variances[order(variances, decreasing = TRUE)])[1:(length(variances)/4)]
vsn_mat <- vsn_mat[topx, ]
NROW(vsn_mat)

# Just plot to see the result of vsn
hist(variances)
dev.off()
boxplot(vsn_mat)

highly_variable_proteins_mofa<-vsn_mat
write.csv(highly_variable_proteins_mofa,'highly_variable_proteins_mofa.csv')
