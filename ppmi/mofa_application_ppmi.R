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
library(ggpubr)

library(dplyr)

outdir_orig='ppmi/plots/'
output_files<- 'ppmi/output/'

source('bladder_cancer/preprocessing.R')
#source('preprocessing.R')
source('ppmi/deseq2_vst_preprocessing_mirnas.R')



### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')





# prerequisites: mass spec preprocessing and desq2 preprocessing
N_FACTORS=10
VISIT='BL'
VISIT='V08'
VISIT='BL'
TISSUE='CSF'; 
TISSUE='Plasma';

TOP_PN=0.90
TOP_GN=0.10# 0.20
TOP_MN=0.50

MIN_COUNT_G=100
MIN_COUNT_M=10

#TISSUE='untargeted'

metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)



p_params<- paste0(VISIT, '_', TISSUE, '_', TOP_PN, '_')
g_params<-paste0(VISIT, '_', TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(VISIT, '_', TOP_MN, '_', MIN_COUNT_M, '_') 
mofa_params<-paste0(N_FACTORS )
#


highly_variable_proteins_outfile = paste0(output_files, p_params , 'highly_variable_proteins_mofa.csv')
fname<-paste0(output_files, 'proteomics_bl.csv')


highly_variable_genes_outfile<-paste0(output_files, 'rnas_',g_params,'_highly_variable_genes_mofa.csv')
highly_variable_mirnas_outfile<-paste0(output_files, 'mirnas_',m_params,'_highly_variable_genes_mofa.csv')



out_params<- paste0( 'p_', p_params, 'g_', g_params, 'm_', m_params, mofa_params )
outdir = paste0(outdir_orig,out_params , '/');outdir

dir.create(outdir, showWarnings = FALSE)
fname<-paste0(output_files, 'proteomics_', VISIT, '_',TISSUE, '.csv')
fname





outdir
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
highly_variable_proteins_mofa<-selectMostVariable(proteomics, 0.15)
dim(highly_variable_proteins_mofa)
#highly_variable_proteins_mofa<-proteomics

###### Option 2: preprocessing from VSN: Load from file 
### 1. Load proteins 

in_file<-highly_variable_proteins_outfile

highly_variable_proteins_mofa<-as.matrix(fread(in_file,header=TRUE), rownames=1)





### Start loading mofa data
proteomics<-as.data.frame(highly_variable_proteins_mofa)
dim(proteomics)


##### Load mirnas + RNAs 
highly_variable_mirnas_mofa<-fread(highly_variable_mirnas_outfile,header=TRUE)
colnames(highly_variable_mirnas_mofa)[1]<-'mirnas'
rownames(highly_variable_mirnas_mofa)<-highly_variable_mirnas_mofa$mirnas

# or input to vst or put as is normalized

miRNA<-as.data.frame(highly_variable_mirnas_mofa[, mirnas:=NULL])
rownames(miRNA)<-rownames(highly_variable_mirnas_mofa)
head(rownames(miRNA))



##### Load RNA seq: 

highly_variable_genes_mofa<-fread(highly_variable_genes_outfile,header=TRUE)
colnames(highly_variable_genes_mofa)[1]<-'rnas'
rownames(highly_variable_genes_mofa)<-highly_variable_genes_mofa$rnas

# or input to vst or put as is normalized

RNA<-as.data.frame(highly_variable_genes_mofa[, rnas:=NULL])
rownames(RNA)<-rownames(highly_variable_genes_mofa)
head(rownames(RNA))

dim(RNA)


##### duplicates 
## histograms to check normal pattern 
library(ggplot2)
miRNAm<-melt(miRNA)
proteomicsm<-melt(proteomics)
par(mfrow=c(1,3))
RNAm<-melt(RNA)

p1<-ggplot(miRNAm, aes(x=value))+ geom_histogram()+ labs(title='mirnas')
ggsave(paste0(outdir, 'data_histograms_mirnas.jpeg' ), width = 10, height=8)

p2<-ggplot(RNAm, aes(x=value))+ geom_histogram()+ labs(title='rnas')
ggsave(paste0(outdir, 'data_histograms_rnas.jpeg' ), width = 10, height=8)

p3<-ggplot(proteomicsm, aes(x=value))+ geom_histogram()+ labs(title='proteins')
ggsave(paste0(outdir, 'data_histograms_proteins.jpeg' ), width = 10, height=8)

dev.off()
dim(highly_variable_proteins_mofa)

dim(proteomics)

##### Filter samples that have all modalities present!

common_samples<-intersect(colnames(miRNA), colnames(proteomics))
common_samples<-intersect(common_samples,combined_bl$PATNO)
common_samples<-intersect(common_samples,colnames(RNA))


common_samples


#### Filter samples that are common in all three
ids<-rownames(miRNA)
miRNA<-as.data.frame(miRNA)
miRNA_filt<-miRNA[which(colnames(miRNA) %in% common_samples)]
miRNA_filt<-miRNA_filt[ ,common_samples]
rownames(miRNA_filt)
#miRNA_filt<-miRNA %>% select(common_samples)


### Select creates the new matrix with the same order 
prot_filt<-proteomics[which(colnames(proteomics) %in% common_samples)]
prot_filt<-prot_filt[, common_samples]


RNA_filt<-RNA[which(colnames(RNA) %in% common_samples)]
RNA_filt<-RNA_filt[, common_samples]

rownames(RNA_filt)
#RNA_filt<-as.data.table(RNA) %>% select(common_samples)





rownames(prot_filt)<-rownames(proteomics)
dim(miRNA_filt)
dim(prot_filt)


## how to analyze and normalize olink? 

mat<-as.matrix(miRNA_filt)


##### Extract metadata with PATNO index






### INPUT TO MOFA
head(unique(rownames(miRNA_filt)))
head(rownames(RNA_filt))
head(rownames(prot_filt))

data = list(proteomics = as.matrix(prot_filt),
            miRNA=as.matrix(miRNA_filt), 
            RNA=as.matrix(RNA_filt) )




###  REMOVE NON ens ids 

order_cols<-colnames(miRNA_filt)
colnames(prot_filt)

NCOL(miRNA_filt)


######
# Add metadata
combined_bl$X<-NULL
match( common_samples, combined_bl$PATNO)
metadata_filt<-combined_bl[match( common_samples, combined_bl$PATNO), ]
NROW(unique(metadata_filt$PATNO))


##### Setup MOFA model 
## model opts 
# SET factor values 
MOFAobject <- create_mofa(data)
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- N_FACTORS
model_opts
MOFAobject <- prepare_mofa(MOFAobject,
                           model_options = model_opts
)


plot_data_overview(MOFAobject)
outdir
ggsave(paste0(outdir, 'data_overview.jpeg'))


##### run the model 
mofa_file<-paste0(outdir,'mofa_ppmi.hdf5')
if (file.exists(mofa_file)){
  pre_trained<-load_model(paste0(outdir,'mofa_ppmi.hdf5'))
  MOFAobject<-pre_trained
  
  
}else {
  MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'))
  
}

###


##### Basic stats

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 4, height=4, dpi=100)


#MOFAobject@training_stats
######
# Check model metadata
metadata_filt$sample<-as.character(metadata_filt$PATNO)
samples_metadata(MOFAobject)<-metadata_filt
#samples_metadata(MOFAobject)
MOFAobject



