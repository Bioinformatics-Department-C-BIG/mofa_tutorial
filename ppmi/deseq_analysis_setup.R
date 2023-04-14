
script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))

print(script_dir)
suppressWarnings(library('R.filesets' ))
suppressWarnings(library(DESeq2))
suppressWarnings(library("SummarizedExperiment"))
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))


### TODO: Add volcano plot for each time point -DONE
### TODO: add heatmap for all tps tpogether -DONE
#source('ppmi/de')

#load
suppressWarnings(library("factoextra"))
suppressWarnings(library("FactoMineR"))
suppressWarnings(library('pheatmap'))
suppressWarnings(library('ggplot2'))

#### Run DE


source(paste0(script_dir, '/../bladder_cancer/preprocessing.R'))
source(paste0(script_dir, '/utils.R'))

VISIT='V08'


process_mirnas<-FALSE


source(paste0(script_dir, '/config.R'))



print(deseq_file)

datalist=loadRDS(deseq_file)
ddsSE=datalist[[1]]
vsd=datalist[[2]]
se_filt=datalist[[3]]
deseq2Results=datalist[[4]]

table(se_filt$COHORT_DEFINITION)

# todo join strings
# TODO: Report the number of samples too! 

des<-gsub(' ', '', paste0(as.character(design(ddsSE))[-1]))

#if  (process_mirnas){
#
#  outdir_s<-paste0(outdir_orig, '/single/', param_str_m_f, des)
#  
#}else{
#  outdir_s<-paste0(outdir_orig, '/single/', param_str_g_f, des)
#  
#}
#outdir_s
# RUN DIFFERENTIAL EXPRESSION ANALYSIS 


### TODO: save the file so we don't have to fit the model each time!! 
dds<-ddsSE



#rm(deseq2Data)



suppressWarnings(dir.create(outdir_s))
pca_files<-paste0(outdir_s, '/PCA/')

#deseq2Data<-loadRDS(paste0(outdir_s, '/deseq_results.RDS'))
#### First obtain the single omics significant RNAs 

