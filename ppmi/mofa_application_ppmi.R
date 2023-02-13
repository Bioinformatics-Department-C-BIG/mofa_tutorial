if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")



#BiocManager::install("MOFA2")
#devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
#browseVignettes("MOFA2")
#BiocManager::install("MOFAdata")

library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggplot2)




source('bladder_cancer/preprocessing.R')
source('ppmi/deseq2_vst_preprocessing_mirnas.R')


combined_bl<-read.csv2('combined_bl.csv', row.names = FALSE)


### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')


# prerequisites: mass spec preprocessing and desq2 preprocessing



outdir='ppmi/plots/'
output_files<- 'ppmi/output/'


csf=1;untargeted=0
fname<-paste0(output_files, 'proteomics_bl.csv')
if (csf){
  outdir=paste0(outdir, '/csf/')
  output_files_prot=paste0(output_files, '/csf/')
  fname<-paste0(output_files_prot, 'proteomics_csf_bl.csv')
  
}  else if(untargeted){
  outdir=paste0(outdir, '/untargeted/')
  output_files_prot=paste0(output_files, '/untargeted/')
  fname<-paste0(output_files_prot, 'untargeted_prot_bl.csv')
  
}


##### Create MOFA object 

##### Handle duplicates 
#proteomics<-proteomics[,!duplicated(colnames(proteomics),fromLast=TRUE)]
output_files

prot_bl_wide<-fread(fname, header=TRUE)
colnames(prot_bl_wide)[1]<-'protein'

proteomics<-as.data.frame(prot_bl_wide, row.names=prot_bl_wide$protein)
rownames(proteomics)<-proteomics$protein
proteomics$protein<-NULL
hist(as.numeric(as.matrix(proteomics)))


# Extract highly variable proteins and genes 
highly_variable_proteins_mofa<-selectMostVariable(proteomics, 0.3)
#highly_variable_proteins_mofa<-proteomics


untargeted=0;csf=1;

###### Option 2 preprocessing from VSN: Load from file 
if (csf){
  in_file<-paste0(output_files_prot, 'highly_variable_proteins_mofa.csv')
}else if(untargeted){
  in_file<-paste0(output_files_prot, 'highly_variable_proteins_mofa.csv')
  
  
}

in_file

highly_variable_proteins_mofa<-as.matrix(fread(in_file,header=TRUE), rownames=1)

proteomics<-highly_variable_proteins_mofa
dim(proteomics)

highly_variable_genes_mofa<-fread(paste0(output_files, 'highly_variable_genes_mofa.csv'),header=TRUE)
colnames(highly_variable_genes_mofa)[1]<-'mirnas'
rownames(highly_variable_genes_mofa)<-highly_variable_genes_mofa$mirnas

# or input to vst or put as is normalized

miRNA<-highly_variable_genes_mofa[, mirnas:=NULL]
rownames(miRNA)
#miRNA<- mirnas_BL

##### duplicates 
## histograms to check normal pattern 
library(ggplot2)
miRNAm<-melt(miRNA)
proteomicsm<-melt(proteomics)
ggplot(miRNAm, aes(x=value))+ geom_histogram()
ggsave(paste0(outdir, 'mirnas_hist.png' ))
ggplot(proteomicsm, aes(x=value))+ geom_histogram()
ggsave(paste0(outdir, 'proteins_hist.png' ))


dim(highly_variable_proteins_mofa)

dim(proteomics)

common_samples<-intersect(colnames(miRNA), colnames(proteomics))
common_samples<-intersect(common_samples,combined_bl$PATNO)

common_samples
which(colnames(miRNA) %in% common_samples)

NROW(unique(common_samples))
miRNA_filt<-miRNA  %>% select(common_samples)
### Select creates the new matrix with the same order 
prot_filt<-as.data.table(proteomics) %>% select(common_samples)
rownames(prot_filt)<-rownames(proteomics)
dim(miRNA_filt)
dim(prot_filt)


## how to analyze and normalize olink? 

mat<-as.matrix(miRNA_filt)


##### Extract metadata with PATNO index






### INPUT TO MOFA
unique(rownames(miRNA_filt))
rownames(miRNA_filt)
rownames(prot_filt)

data = list(proteomics = as.matrix(prot_filt, rownames=rownames(prot_filt)),
            miRNA=as.matrix(miRNA_filt, rownames=rownames(miRNA_filt)))
ss



###  REMOVE NON ens ids 

order_cols<-colnames(miRNA_filt)
colnames(prot_filt)

NCOL(miRNA_filt)

hist(data$miRNA)
hist(data$proteomics)
dev.off()
plot_data_overview(MOFAobject)

######
# Add metadata
combined_bl$X<-NULL
match( common_samples, combined_bl$PATNO)
metadata_filt<-combined_bl[match( common_samples, combined_bl$PATNO), ]

metadata_filt$PATNO
NROW(unique(metadata_filt$PATNO))


#### model opts 
# SET factor values 
MOFAobject <- create_mofa(data)

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 8
model_opts
MOFAobject <- prepare_mofa(MOFAobject,
                           model_options = model_opts
)


##### run the model 

MOFAobject <- run_mofa(MOFAobject, outfile = paste0(output,'mofa_ppmi.hdf5'))

#pre_trained<-load_model(paste0(output, 'mofa_ppmi.hdf5'))
#MOFAobject<-pre_trained
###


##### Basic stats
plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)




######
# Check model metadata
metadata_filt$sample<-as.character(metadata_filt$PATNO)
samples_metadata(MOFAobject)<-metadata_filt
samples_metadata(MOFAobject)


