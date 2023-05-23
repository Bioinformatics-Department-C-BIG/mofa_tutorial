#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")


isRStudio <- Sys.getenv("RSTUDIO") == "1"
if (isRStudio){
  script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
}else{
  script_dir<- "D:/DATADRIVE/Efi Athieniti/Documents/git/mofa/ppmi"
  
}

source(paste0(script_dir,'/setup_os.R'))


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
library(ggplot2)


library('MultiAssayExperiment')


outdir_orig=paste0(data_dir,'ppmi/plots/')
output_files<- paste0(data_dir,'ppmi/output/')

source(paste0(script_dir,'/../bladder_cancer/preprocessing.R'))

#source('preprocessing.R')
#source('ppmi/deseq2_vst_preprocessing_mirnas.R')



### preprocessing 
#1. dsb proteomics 2. RNseq: just log
#install.packages('psych')





# prerequisites: mass spec preprocessing and desq2 preprocessing

# TODO: move all to config file 
split=FALSE
run_rna_mirna=FALSE
### if we are using all modalities we might need to change TOP_GN
FULL_SET=TRUE
VISIT_COMPARE='BL'
# cohort 1 =prodromal 

if (split){
  N_FACTORS=8
}
VISIT=c('V08');
run_vsn=TRUE

## tissue is set in the config


NORMALIZED=TRUE;
use_signif=FALSE
process_mirnas=FALSE

source(paste0(script_dir, '/config.R'))
source(paste0(script_dir, '/mofa_config.R'))

NORMALIZED
TOP_GN; TOP_MN
metadata_output<-paste0(output_files, 'combined.csv')
combined_all_original<-read.csv2(metadata_output)

metadata_output<-paste0(output_files, 'combined_log.csv')
combined<-read.csv2(metadata_output)
dim(combined)

combined_bl<-combined
which(is.na(combined_bl$AGE))
combined_bl$AGE
scale_views=TRUE
run_mofa_complete<-FALSE

## MORE samples, more factors ! 
# TODO: need a better way to decide how to run this 
if (run_mofa_complete){
  N_FACTORS=9
  
}else{
  N_FACTORS=12
  
}


#combined$Outcome
## VISIT_S to allow this to be more than one visits at once!! 

#p_params<- paste0(VISIT_S, '_',TISSUE, '_', TOP_PN, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)



mofa_params<-paste0(N_FACTORS,'_sig_',  use_signif,'complete', run_mofa_complete )

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
highly_variable_genes_outfile
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


mofa_params;g_params;p_params;p_params_out
p_params
out_params<- paste0( 'p_', p_params, 'g_', g_params, 'm_', m_params, mofa_params, '_coh_', sel_coh_s,'_', VISIT_S, '_', scale_views[1])
highly_variable_proteins_outfile<-paste0(output_files, p_params , '_highly_variable_proteins_mofa.csv')
out_params

outdir = paste0(outdir_orig,out_params, '_split_', split , '/');outdir
dir.create(outdir, showWarnings = FALSE)
outdir
fname<-paste0(output_files, 'proteomics_',TISSUE, '.csv')




###### Option 2: preprocessing from VSN: Load from file 
### 1. Load proteins 

in_file<-highly_variable_proteins_outfile

highly_variable_proteins_mofa<-as.matrix(fread(in_file,header=TRUE), rownames=1)

### Start loading mofa data
proteomics<-as.data.frame(highly_variable_proteins_mofa)
dim(highly_variable_proteins_mofa)

test2<-proteomics['GLO1',]
##### Load mirnas + RNAs 
### we use data.table because there are duplicate samples? 
### problem with saving of rownmaes 
highly_variable_mirnas_mofa<-fread(highly_variable_mirnas_outfile,header=TRUE)
colnames(highly_variable_mirnas_mofa)[1]<-'mirnas'
rownames(highly_variable_mirnas_mofa)<-highly_variable_mirnas_mofa$mirnas


highly_variable_mirnas_outfile
# EITHER input to vst or put as is normalized
miRNA<-as.data.frame(highly_variable_mirnas_mofa[, mirnas:=NULL])
rownames(miRNA)<-rownames(highly_variable_mirnas_mofa)
head(rownames(miRNA));
head(colnames(miRNA))



##### Load RNA seq: 

highly_variable_genes_mofa<-fread(highly_variable_genes_outfile,header=TRUE)
colnames(highly_variable_genes_mofa)[1]<-'rnas'
rownames(highly_variable_genes_mofa)<-highly_variable_genes_mofa$rnas
dim(highly_variable_genes_mofa)
# or input to vst or put as is normalized

RNA<-as.data.frame(highly_variable_genes_mofa[, rnas:=NULL])
rownames(RNA)<-rownames(highly_variable_genes_mofa)
head(rownames(RNA)); head(colnames(RNA))



create_hist<-function(df, name){
  
  dfm<-melt(df)
  
  p1<-ggplot(dfm, aes(x=value))+ geom_histogram()+ labs(title='mirnas')
  ggsave(paste0(outdir, 'data_histograms',name,  '.jpeg' ), width = 10, height=8)
}
## histograms to check normal pattern 
create_hist(RNA, 'RNA')
create_hist(miRNA, 'miRNA')




############################# Preprocessing ############################
########################################################################
### Create the summarized experiment 
########
########


#data = list(
#  miRNA=as.matrix(miRNA_filt), 
#  RNA=as.matrix(RNA_filt) )

#### just trying a multi assay here to help with filtering.. 
#head(colnames(prot_filt));head(colnames(miRNA_filt)); colnames(RNA_filt)

### might need to filter by what is common with meta


data_full<-list(miRNA=as.matrix(miRNA), 
                RNA=as.matrix(RNA),
                proteomics=as.matrix(proteomics))


assay_full=c(rep('RNA', length(RNA)),
             rep('miRNA', length(miRNA)),
             rep('proteomics', length(proteomics)))

colname = c(colnames(RNA), colnames(miRNA), colnames(proteomics))
primary=colname
sample_map=DataFrame(assay=assay_full, primary=primary, colname=colname)
common_samples_in_assays=unique(colname)
### TODO: is it a problem for duplicates when i make patno_event_id the key column? 
### Note: HERE WE lost duplicate metadata ie. double clinical measures for one patient

metadata_filt<-combined_bl[match(common_samples_in_assays, combined_bl$PATNO_EVENT_ID),]
metadata_filt$primary<-metadata_filt$PATNO_EVENT_ID

rownames(metadata_filt)=metadata_filt$PATNO_EVENT_ID


mofa_multi<-MultiAssayExperiment(experiments=data_full,
                                 colData = metadata_filt, 
                                 sampleMap=sample_map)


mofa_multi_complete_all<-mofa_multi[,complete.cases(mofa_multi)]
mofa_multi_complete_all<-mofa_multi[,complete.cases(mofa_multi)]

mofa_multi_rna_mir<-subsetByAssay(mofa_multi, c('RNA', 'miRNA'))
mofa_multi_rna_mir_complete<-mofa_multi_rna_mir[,complete.cases(mofa_multi_rna_mir)]

library('UpSetR')
upsetSamples(mofa_multi)
#mofa_multi_V04=mofa_multi[,mofa_multi$EVENT_ID %in% VISIT]


###  REMOVE NON ens ids 

nsamples<-dim(colData(mofa_multi_complete_all))[1]
nsamples


### Split the data
if (split){
  seed_tr_test=150
  set.seed(seed_tr_test)
  train_ind<-sample(nsamples, nsamples*0.7)
  mofa_multi_complete_train = mofa_multi_complete_all[,train_ind]
  mofa_multi_complete_test = mofa_multi_complete_all[,-train_ind]
  mofa_multi_complete=mofa_multi_complete_train
}else{
  mofa_multi_complete=mofa_multi_complete_all
  mofa_multi_complete=mofa_multi_complete_all
  
}
#dim(colData(mofa_multi_complete_train))[1]
prot_to_impute<-assays(mofa_multi_complete)$proteomics


###################### RUN MOFA #########################
##### Setup MOFA model 
## model opts 
# SET factor values 



#N_FACTORS=8
### separate visits 
outdir

if (run_mofa_complete){
  MOFAobject <- create_mofa(mofa_multi_complete)
  
}else{
  MOFAobject <- create_mofa(mofa_multi)
  
}

mofa_multi

if (length(VISIT)>1){
  MOFAobject <- create_mofa(mofa_multi_complete, groups= mofa_multi_complete$EVENT_ID)
  
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




outdir
ggsave(paste0(outdir, 'data_overview.jpeg'))

#### TODO FIX THE DATAFRAME 
#outdir = paste0(outdir_orig,out_params , '_', VISIT, '/');
outdir
dir.create(outdir, showWarnings = FALSE)
##### run the model 

#MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'), use_basilisk = TRUE)


mofa_file<-paste0(outdir,'mofa_ppmi.hdf5')
if (file.exists(mofa_file)){
  pre_trained<-load_model(paste0(outdir,'mofa_ppmi.hdf5'))
  MOFAobject<-pre_trained
  
  
}else {
  MOFAobject <- run_mofa(MOFAobject, outfile = paste0(outdir,'mofa_ppmi.hdf5'), use_basilisk = TRUE)
}
##### Basic stats

plot_variance_explained(MOFAobject, max_r2=20)
ggsave(paste0(outdir, 'variance_explained_total','.png'), width = 7, height=4, dpi=100)


samples_metadata(MOFAobject)$Outcome





