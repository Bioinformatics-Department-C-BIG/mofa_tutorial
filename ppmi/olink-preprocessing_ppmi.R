
#BiocManager::install('DEP')
## TODO: change all scripts to be agnostic of visit until mofa
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library(ggplot2)
library("vsn")
library("DEP")
library("data.table")
library("SummarizedExperiment")
script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)


if (!require("pacman")) install.packages("pacman")
#pacman::p_load(dplyr,tidyr,DESeq2,edgeR,limma,ComplexHeatmap,EnhancedVolcano,tibble,fgsea,stringr,org.Hs.eg.db)
source(paste0(script_dir,'/../bladder_cancer/preprocessing.R'))

output_1=paste0('ppmi/plots/')
output_files<-paste0('ppmi/output/')


TOP_PN<-0.7

param_str<-paste0(TOP_PN)



TISSUE='CSF'

TISSUE='CSF'


VISIT='BL'

VISIT='BL'
VISIT='BL'
TISSUE='Plasma'
NORMALIZED=TRUE

#### read in proteomics 
p_params_in<- paste0( TISSUE, '_', NORMALIZED)
p_params_out<- paste0(TISSUE, '_', TOP_PN, '_', NORMALIZED)


if (NORMALIZED){
  in_file_original<-paste0(output_files, 'proteomics_', p_params_in,  '_no_log.csv')
  
  
}else{
  in_file_original<-paste0(output_files, 'proteomics_', p_params_in, '.csv')


}

highly_variable_proteins_outfile<-paste0(output_files, p_params_out , '_highly_variable_proteins_mofa.csv')

# Read in 
prot_bl_wide_unlog<-as.matrix(fread(in_file_original, header=TRUE), rownames=1)
proteomics<-prot_bl_wide_unlog
  



# Remove rows with 90% NA 
df<-proteomics
proteomics <- df[rowSums(is.na(df)) < round(0.2*ncol(df)), ]
dim(df); dim(proteomics)


### filter here before editing more 
## filter out rows with very low min count
df<-proteomics; dim(df)
min.count= quantile(df, na.rm = TRUE, 0.1)
min.count= min(df, na.rm = TRUE)
dim(df)
keep <- rowSums(df>min.count, na.rm = TRUE) >= round(0.9*ncol(df))
length(which(keep))
proteomics<-proteomics[keep,]
dim(df); dim(proteomics)



dim(proteomics)
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

assays(se)[[1]][1]
data_columns
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
head(rownames(vsn_mat))
head(colnames(vsn_mat))
hist(vsn_mat)

boxplot(vsn_mat[,1:30])
dim(normalized_data)

# Select the top most variable proteins
## TODO: fix the bug in selectMostVariable

highly_variable_proteins_mofa2<-selectMostVariable(vsn_mat, TOP_PN)


variances<-rowVars(vsn_mat, na.rm = TRUE)
maxn<-round(length(variances)*TOP_PN)
to_sel<-rownames(vsn_mat)[order(variances, decreasing=TRUE)][1:maxn]
rownames(vsn_mat[to_sel,])
var(vsn_mat['NPPB',], na.rm=TRUE)
highly_variable_proteins_mofa=vsn_mat[to_sel,]

highly_variable_proteins_mofa2

dim(highly_variable_proteins_mofa)
rownames(highly_variable_proteins_mofa)
# Just plot to see the result of vsn
boxplot(highly_variable_proteins_mofa[,1:70])

colnames(highly_variable_proteins_mofa)<-sample

write.csv(highly_variable_proteins_mofa,highly_variable_proteins_outfile)
dim(highly_variable_proteins_mofa)

hist(highly_variable_proteins_mofa)





