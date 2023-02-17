
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
library('data.table')

library('SummarizedExperiment')
library('ggplot2')


if (!require("pacman")) install.packages("pacman")
#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source('bladder_cancer/preprocessing.R')

untargeted=1;csf=0
outdir='ppmi/plots/'
output_files<-'ppmi/output/'
if (csf){
  outdir=paste0(outdir, '/csf/')
}  else if(untargeted){
  outdir=paste0(outdir, '/untargeted/')
  
} else if(plasma){
  outdir=paste0(outdir, '/plasma/')
}


#### read in proteomics 
in_file<-paste0(output_files, 'untargeted_prot_bl.csv')
prot_bl_wide<-as.matrix(fread(in_file, header=TRUE), rownames=1)
proteomics<-prot_bl_wide






##### Remove high NA







filter_nas<-function(data){
    data<-as.data.frame(data)
    find_na<-apply(is.na(data),1,sum)
    thresh<-round(dim(data)[2]/6)
    # only select ids that have less than x in na
    ids<-which(find_na<thresh)
    data_filt<-data[ids,]
    return(data_filt)
}

proteomics<-filter_nas(proteomics)
data<-proteomics

data$name<-c(rownames(data))
data$ID<-data$name

data$name
data$ID
data$ID
data_columns=seq(1:dim(proteomics)[2])
data_columns

sample<-colnames(proteomics)
sample






##### Make SummarizedExperiment


###todo 
exp_design = data.frame(label=sample,condition=sample, replicate=rep(1, dim(proteomics)[2]))
exp_design




se <- make_se(data, data_columns,exp_design)



###
#plot_missval(se)


is.nan(as.matrix(data))
boxplot(log(data[1:50]))

ggsave(paste0(outdir,'boxplot.png'))

interm<-as.matrix(assays(se)[[1]])

is.nan(as.matrix(interm))


# Filter and normalize
normalized_data<-normalize_vsn(se)


meanSdPlot(normalized_data)
ggsave(paste0(outdir,'meansd.png' ))
# Check plot after vsn
vsn_mat<-assays(normalized_data)[[1]]




hist(as.numeric(as.matrix(vsn_mat)))


boxplot(vsn_mat[,1:30])
dim(normalized_data)

# Select the top most variable proteins
highly_variable_proteins_mofa<-selectMostVariable(vsn_mat, 0.25)
dim(highly_variable_proteins_mofa)
rownames(highly_variable_proteins_mofa)
# Just plot to see the result of vsn
#boxplot(highly_variable_proteins_mofa)

colnames(highly_variable_proteins_mofa)<-sample

write.csv(highly_variable_proteins_mofa,paste0(output_files,'/untargeted/highly_variable_proteins_mofa.csv'))


