if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir,'/setup_os.R'))

#BiocManager::install('MultiAssayExperiment')
BiocManager::install("MOFA2")
#devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"), force = TRUE)
#browseVignettes("MOFA2")
#BiocManager::install("MOFAdata")'

# TODO: update to receive the whole matrices and start filtering here depending on the test we want to do
# SCENARIOS: 
# select cohort: 1,2,3,4: PD, Prodromal, , Healthy Control
# select visit: ALL, V02, V04, V06, V08 
# question: should we preprocces each cohort separately?


library(MOFA2)
library(data.table)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(dplyr)

library('MultiAssayExperiment')


outdir_orig=paste0(data_dir,'ppmi/plots/')
output_files<- paste0(data_dir,'ppmi/output/')

source(paste0(script_dir,'/../bladder_cancer/preprocessing.R'))


source(paste0(script_dir,'/config.R'))
#source('preprocessing.R')
#source('ppmi/deseq2_vst_preprocessing_mirnas.R')



### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')





# prerequisites: mass spec preprocessing and desq2 preprocessing

# TODO: move all to config file 

TOP_PN=0.70

FULL_SET=TRUE

NA_PERCENT=0.8


VISIT_COMPARE='BL'
TOP_PN=0.90
# cohort 1 =prodromal 
sel_coh <- c(1,4)

VISIT='BL'


VISIT=c('V08')

VISIT=c('BL')
VISIT=c('V04')
#VISIT=c('BL')
VISIT=c('BL')

VISIT=c('V06');
sel_coh <- c(2)


N_FACTORS=8

VISIT=c('V08');





TISSUE='CSF'; 
run_vsn=TRUE
TISSUE='Plasma';


NORMALIZED=TRUE;
use_signif=FALSE

source(paste0(script_dir, '/config.R'))

metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)
combined_bl<-combined

scale_views=TRUE

combined$Outcome
## VISIT_S to allow this to be more than one visits at once!! 

p_params<- paste0(VISIT_S, '_',TISSUE, '_', TOP_PN, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)

mofa_params<-paste0(N_FACTORS,'_sig_',  use_signif )
#param_str_g<-paste0('rnas_', g_params, sel_coh_s, '_'  )
#


highly_variable_proteins_outfile = paste0(output_files, p_params , '_highly_variable_proteins_mofa.csv')
highly_variable_genes_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
highly_variable_mirnas_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')

if (use_signif){
  highly_variable_genes_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa_signif.csv')
  highly_variable_mirnas_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa_signif.csv')
  
}

highly_variable_mirnas_outfile
highly_variable_proteins_outfile

### TODO: INPUT vsn files and filter here instead of rerunnign!!1 
#param_str_m<-paste0('mirnas_', m_params ,sel_coh_s, '_')
#vsn_file_m<-paste0(output_files, 'mirnas_', param_str_m, '_vsn.csv')
#vsn_file_g<-paste0(output_files, 'rnas_', param_str_g, '_vsn.csv')

###
#vsn_file_m_df<-read.csv(vsn_file_m,header=TRUE)
#duplicated(colnames(vsn_file_m_df))
#vsn_file_m_df[,1]
#rownames(vsn_file_m_df)<-vsn_file_m_df[,1]



out_params<- paste0( 'p_', p_params, 'g_', g_params, 'm_', m_params, mofa_params, '_coh_', sel_coh_s,'_', VISIT_S, '_', scale_views[1])
outdir = paste0(outdir_orig,out_params , '/');outdir
dir.create(outdir, showWarnings = FALSE)

fname<-paste0(output_files, 'proteomics_',TISSUE, '.csv')
fname



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


highly_variable_mirnas_outfile
# EITHER input to vst or put as is normalized
miRNA<-as.data.frame(highly_variable_mirnas_mofa[, mirnas:=NULL])
rownames(miRNA)<-rownames(highly_variable_mirnas_mofa)
head(rownames(miRNA))



##### Load RNA seq: 

highly_variable_genes_mofa<-fread(highly_variable_genes_outfile,header=TRUE)
colnames(highly_variable_genes_mofa)[1]<-'rnas'
rownames(highly_variable_genes_mofa)<-highly_variable_genes_mofa$rnas
dim(highly_variable_genes_mofa)
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
colnames(proteomics)
##### Filter samples that have all modalities present!


common_samples<-intersect(colnames(miRNA), colnames(proteomics)); common_samples

## do not add proteomics
common_samples<-intersect(colnames(miRNA), colnames(miRNA)); common_samples
common_samples<-intersect(common_samples,colnames(RNA)) ; common_samples
common_samples<-intersect(common_samples,combined_bl$PATNO_EVENT_ID); common_samples

######
# Add metadata


metadata_filt<-combined_bl[match(common_samples, combined_bl$PATNO_EVENT_ID),]
only_pd<-metadata_filt$PATNO_EVENT_ID[which(metadata_filt$COHORT %in% sel_coh)]


common_samples<-only_pd;common_samples
metadata_filt<-metadata_filt[match(only_pd, metadata_filt$PATNO_EVENT_ID),]
# Rewrite to add only pd
#metadata_filt[c('COHORT', 'COHORT_DEFINITION')]

#metadata_filt$PATNO

write.csv(common_samples,paste0(output_files,out_params, '_common_samples.txt'), 
          row.names = FALSE, quote=FALSE )

head(common_samples)


#### Filter samples that are common in all three
ids<-rownames(miRNA)
miRNA<-as.data.frame(miRNA); dim(miRNA)
miRNA_filt<-miRNA[,match(common_samples, colnames(miRNA)) ]
dim(miRNA_filt)
#miRNA_filt<-miRNA_filt[ ,common_samples]


### Select creates the new matrix with the same order 
prot_filt<-proteomics[,match(common_samples, colnames(proteomics))]
dim(prot_filt)

# use match to get column too
RNA_filt<-RNA[,match(common_samples, colnames(RNA)) ]
dim(prot_filt)


head(rownames(RNA_filt))
head(rownames(miRNA_filt))
head(rownames(prot_filt))

#RNA_filt<-as.data.table(RNA) %>% select(common_samples)


rownames(prot_filt)<-rownames(proteomics)
dim(miRNA_filt)
dim(prot_filt)
#colnames(miRNA_filt)
#colnames(RNA_filt)


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

data = list(
            miRNA=as.matrix(miRNA_filt), 
            RNA=as.matrix(RNA_filt) )

#### just trying a multi assay here to help with filtering.. 

head(colnames(prot_filt));head(colnames(miRNA_filt)); colnames(RNA_filt)

assay=c(rep('proteomics', length(prot_filt)),
        rep('miRNA', length(miRNA_filt)),
        rep('RNA', length(RNA_filt)))

assay=c(rep('miRNA', length(miRNA_filt)),
        rep('RNA', length(RNA_filt)))


primary=metadata_filt$PATNO_EVENT_ID
colname=metadata_filt$PATNO_EVENT_ID
sample_map=DataFrame(assay=assay, primary=primary, colname=colname)
rownames(metadata_filt)=metadata_filt$PATNO_EVENT_ID

mofa_multi<-MultiAssayExperiment(experiments=data,
                     colData = metadata_filt, 
                     sampleMap=sample_map)

complete.cases(metadata_filt$EVENT_ID)
#install.packages('UpSetR')
library('UpSetR')
upsetSamples(mofa_multi)
mofa_multi_V04=mofa_multi[,mofa_multi$EVENT_ID %in% VISIT]

mofa_multi_V04

###  REMOVE NON ens ids 

order_cols<-colnames(miRNA_filt)
head(cbind(colnames(prot_filt),colnames(RNA_filt), colnames(miRNA_filt)  ))

NCOL(miRNA_filt)

colData(mofa_multi_V04)

##### Setup MOFA model 
## model opts 
# SET factor values 



#N_FACTORS=8
### separate visits 
outdir
MOFAobject <- create_mofa(mofa_multi_V04)

if (length(VISIT)>1){
  MOFAobject <- create_mofa(mofa_multi_V04, groups= mofa_multi_V04$EVENT_ID)
  
}

model_opts <- get_default_model_options(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
model_opts$num_factors <- N_FACTORS
data_opts
data_opts$scale_views=scale_views
MOFAobject <- prepare_mofa(MOFAobject,
                           model_options = model_opts,
                           data_options = data_opts
)




plot_data_overview(MOFAobject)
outdir
ggsave(paste0(outdir, 'data_overview.jpeg'))

#### TODO FIX THE DATAFRAME 
#outdir = paste0(outdir_orig,out_params , '_', VISIT, '/');
outdir
dir.create(outdir, showWarnings = FALSE)
##### run the model 

#MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'))



#if (file.exists(mofa_file)){
# pre_trained<-load_model(paste0(outdir,'mofa_ppmi.hdf5'))
# MOFAobject<-pre_trained
 
 
#}else {
  MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'), use_basilisk = TRUE)
#}
##### Basic stats

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 7, height=4, dpi=100)


samples_metadata(MOFAobject)$Outcome






