if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



#BiocManager::install("MOFA2")
#devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
#browseVignettes("MOFA2")
#BiocManager::install("MOFAdata")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggplot2)


outdir='ppmi/plots/'
output_files<- 'ppmi/output/'

source('bladder_cancer/preprocessing.R')
source('ppmi/deseq2_vst_preprocessing_mirnas.R')

### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')


# prerequisites: mass spec preprocessing and desq2 preprocessing



##### Create MOFA object 

##### Handle duplicates 
#proteomics<-proteomics[,!duplicated(colnames(proteomics),fromLast=TRUE)]


prot_bl_wide<-fread(paste0(output_files, 'proteomics_bl.csv'), header=TRUE)
colnames(prot_bl_wide)[1]<-'protein'

proteomics<-as.data.frame(prot_bl_wide, row.names=prot_bl_wide$protein)
rownames(proteomics)<-proteomics$protein
proteomics$protein<-NULL
hist(as.numeric(as.matrix(proteomics)))
highly_variable_proteins_mofa<-selectMostVariable(proteomics, 0.5)

proteomics<-highly_variable_proteins_mofa


highly_variable_genes_mofa<-fread(paste0(output_files, 'highly_variable_genes_mofa.csv'),header=TRUE)
colnames(highly_variable_genes_mofa)[1]<-'mirnas'
rownames(highly_variable_genes_mofa)<-highly_variable_genes_mofa$mirnas

# or input to vst or put as is normalized

miRNA<-highly_variable_genes_mofa[, mirnas:=NULL]
rownames(miRNA)
#miRNA<- mirnas_BL

##### dups 
library(ggplot2)
miRNAm<-melt(miRNA)
proteomicsm<-melt(proteomics)
ggplot(miRNAm, aes(x=value))+ geom_histogram()
ggsave(paste0(output, 'mirnas_hist.png' ))
ggplot(proteomicsm, aes(x=value))+ geom_histogram()
ggsave(paste0(output, 'proteins_hist.png' ))


dim(proteomics)

common_samples<-intersect(colnames(miRNA), colnames(proteomics))

common_samples
which(colnames(miRNA) %in% common_samples)

miRNA_filt<-miRNA  %>% select(common_samples)
dim(miRNA_filt)
dim(prot_filt)

### Select creates the new matrix with the same order 
prot_filt<-proteomics %>% select(common_samples)

proteomics

## how to analyze and normalize olink? 

mat<-as.matrix(miRNA_filt)


##### Extract metadata with PATNO index

motor_assess_BL_filt<-motor_assess_BL[motor_assess_BL$PATNO %in% common_samples,]





### INPUT TO MOFA
rownames(miRNA_filt)
rownames(prot_filt)

data = list(miRNA=as.matrix(miRNA_filt), 
            proteomics = as.matrix(prot_filt))
#rownames(data$proteomics)<-paste(rownames(data$proteomics),"1",sep="_")



###  REMOVE NON ens ids 

order_cols<-colnames(miRNA_filt)
colnames(prot_filt)

NCOL(miRNA_filt)

hist(data$miRNA)
hist(data$proteomics)


data$proteomics

MOFAobject <- create_mofa(data)
MOFAobject



plot_data_overview(MOFAobject)

######
#Add metadata

metadata_filt<-motor_assess_BL[match( common_samples, motor_assess_BL$PATNO), ]
metadata_filt2<-non_motor_BL[match( common_samples, non_motor_BL$PATNO), ]

metadata_filt<-merge(metadata_filt, metadata_filt2,by='PATNO' )
NROW(unique(metadata_filt$PATNO))


#### model opts 

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8





model_opts
MOFAobject <- prepare_mofa(MOFAobject,
                           model_options = model_opts
)






MOFAobject <- run_mofa(MOFAobject, outfile = paste0(output,'mofa_ppmi.hdf5'))


#pre_trained<-load_model(paste0(output, 'mofa_ppmi.hdf5'))
#MOFAobject<-pre_trained
###

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


MOFAobject@samples_metadata
###
get_weights(MOFAobject)




######


#metadata
metadata_filt$sample<-as.character(metadata_filt$PATNO)


samples_metadata(MOFAobject)<-metadata_filt





