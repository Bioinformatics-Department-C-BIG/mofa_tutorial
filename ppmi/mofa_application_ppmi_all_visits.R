#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")




source(paste0('ppmi/setup_os.R'))

script_dir
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

source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))

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
run_mofa_complete<-FALSE

source(paste0(script_dir, '/ppmi/config.R'))

source(paste0(script_dir, '/ppmi/mofa_config.R'))

# metada source 
metadata_output<-paste0(output_files, 'combined.csv')
combined_all_original<-read.csv2(metadata_output)
metadata_output<-paste0(output_files, 'combined_log.csv')
combined<-read.csv2(metadata_output)

combined_bl<-combined
which(is.na(combined_bl$AGE))

scale_views=TRUE

## MORE samples, more factors ! 
# TODO: need a better way to decide how to run this 
if (run_mofa_complete){
  N_FACTORS=9
}else{
  N_FACTORS=12
  N_FACTORS=15
 # N_FACTORS=20
  

}





#combined$Outcome
## VISIT_S to allow this to be more than one visits at once!! 

#p_params<- paste0(VISIT_S, '_',TISSUE, '_', TOP_PN, '_', substr(NORMALIZED,1,1), '_', sel_coh_s,'vsn_', substr(run_vsn,1,1), 'NA_', NA_PERCENT)
mofa_params<-paste0(N_FACTORS,'_sig_',  use_signif,'complete', run_mofa_complete )
out_params<- paste0( 'p_', p_params, 'g_', g_params, 'm_', m_params, mofa_params, '_coh_', sel_coh_s,'_', VISIT_S, '_', scale_views[1])

### get mofa parameters 
# 1. output directory
# 2. 
outdir = paste0(outdir_orig,out_params, '_split_', split , '/');outdir
dir.create(outdir, showWarnings = FALSE)
prepare_multi_data<-function(p_params, param_str_g, param_str_m, mofa_params){
  

###### Option 2: preprocessing from VSN: Load from file 
### 1. Load proteins 


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

data_full<-prepare_multi_data(p_params, param_str_g, param_str_m, mofa_params)
  
assay_full=c(rep('RNA', length(RNA)),
             rep('miRNA', length(miRNA)),
             rep('proteomics', length(proteomics)))

colname = c(colnames(RNA), colnames(miRNA), colnames(proteomics))
primary=colname
sample_map=DataFrame(assay=assay_full, primary=primary, colname=colname)
common_samples_in_assays=unique(colname)
### TODO: is it a problem for duplicates when i make PATNO_EVENT_ID the key column? 
### Note: HERE WE lost duplicate metadata ie. double clinical measures for one patient

metadata_filt<-combined_bl[match(common_samples_in_assays, combined_bl$PATNO_EVENT_ID),]
metadata_filt$primary<-metadata_filt$PATNO_EVENT_ID

rownames(metadata_filt)=metadata_filt$PATNO_EVENT_ID


mofa_multi<-MultiAssayExperiment(experiments=data_full,
                                 colData = metadata_filt, 
                                 sampleMap=sample_map)


mofa_multi_complete_all<-mofa_multi[,complete.cases(mofa_multi)]

mofa_multi_rna_mir<-subsetByAssay(mofa_multi, c('RNA', 'miRNA'))
mofa_multi_rna_mir_complete<-mofa_multi_rna_mir[,complete.cases(mofa_multi_rna_mir)]

#library('UpSetR')
#p2<-upsetSamples(mofa_multi)  

#ggsave(paste0(outdir, 'intersection.pdf'), plot=p2)
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

#ggsave(paste0(outdir, 'data_overview.jpeg'))

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

mofa_params<-paste0(N_FACTORS,'_sig_',  use_signif,'complete', run_mofa_complete )

cors_both<-get_correlations_with_coh(MOFAobject)
cors_t<-paste(round(cors_pearson[,'CONCOHORT'][sel_factors], digits=2), collapse=', ')

}

df_stats=  c( TOP_PN, TOP_GN, MIN_COUNT_G, TOP_MN, MIN_COUNT_M, mofa_params, sel_coh_s,VISIT_S,  scale_views[1],  use_signif,
              run_mofa_complete, N_FACTORS,cors_t  )


write.table(t(df_stats), paste0(outdir_orig,'all_stats.csv'), append=TRUE,sep=',', col.names = FALSE)





