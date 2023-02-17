
#BiocManager::install('DEP')

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
library("data.table")
library("SummarizedExperiment")


if (!require("pacman")) install.packages("pacman")
#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source('bladder_cancer/preprocessing.R')


output_1='ppmi/plots/'
output_files<-'ppmi/output/'


top_n<-0.5

param_str<-paste0(top_n)

csf=0;plasma=1
#### read in proteomics 
if (csf){
  output_files_prot<-paste0(output_files, '/csf/')
  in_file<-paste0(output_files_prot, 'proteomics_csf_bl_no_log.csv')
}else if (plasma){
  output_files_prot<-paste0(output_files, '/plasma/')
  in_file<-paste0(output_files, 'proteomics_bl_no_log.csv')
  
}

prot_bl_wide_unlog<-as.matrix(fread(in_file, header=TRUE), rownames=1)
proteomics<-prot_bl_wide_unlog

data<-proteomics
data<-as.data.frame(data)

data$name<-c(rownames(data))
data$ID<-data$name

data$name
data$ID
data$ID# Make SummarizedExperiment
data_columns=seq(1:dim(proteomics)[2])
data_columns

sample<-colnames(proteomics)
sample
exp_design = data.frame(label=sample,condition=sample, replicate=rep(1, dim(proteomics)[2]))
exp_design

se <- make_se(data, data_columns,exp_design)

is.nan(as.matrix(data))
boxplot(log(data[1:16]))


interm<-as.matrix(assays(se)[[1]])

is.nan(as.matrix(interm))


# Filter and normalize
normalized_data<-normalize_vsn(se)

meanSdPlot(normalized_data)
ggsave(paste0(outdir,'meansd.png' ))
# Check plot after vsn
vsn_mat<-assays(normalized_data)[[1]]



#### new function
#data$name<-NULL
#data$ID<-NULL



#d1<-is.nan(data[[1]])
#normalized_data<-justvsn(as.matrix(data))



#vsn_mat<-normalized_data


hist(vsn_mat)

boxplot(vsn_mat[,1:30])
dim(normalized_data)

# Select the top most variable proteins
highly_variable_proteins_mofa<-selectMostVariable(vsn_mat, top_n)
dim(highly_variable_proteins_mofa)
rownames(highly_variable_proteins_mofa)
# Just plot to see the result of vsn
boxplot(highly_variable_proteins_mofa[,1:70])

colnames(highly_variable_proteins_mofa)<-sample

write.csv(highly_variable_proteins_mofa,paste0(output_files_prot,'highly_variable_proteins_mofa.csv'))
dim(highly_variable_proteins_mofa)

