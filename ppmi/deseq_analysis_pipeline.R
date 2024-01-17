

source(paste0('ppmi/setup_os.R'))
#source(paste0('ppmi/mofa_application_ppmi_all_visits.R'))

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))

### disconnect from mofa and other scripts 

VISIT=c('BL','V04', 'V06',  'V08');cell_corr = FALSE
VISIT=c('BL','V04', 'V06',  'V08');cell_corr = FALSE

# Pipeline steps 
# 1. run deseq  - with and without correction for cell types
# 2. volcanot plots 
# 3. PCA batch corrected with vsd or lognorm top genes - colour PCA with cohort/updrs 

source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))
source(paste0(script_dir,'ppmi/plotting_utils.R'))
source(paste0(script_dir,'ppmi/utils.R'))




library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)
library(org.Hs.eg.db)
#BiocManager::install('org.Hs.eg.db')



### FIRST LOAD required files  ####
process_mirnas = TRUE;
source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se_mirs=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
se_mirs_norm=load_se_all_visits(input_file = input_file_mirs, combined=combined_bl_log); 

head(assay(se_mirs))
head(assay(se_mirs_norm))

process_mirnas=FALSE
source(paste0(script_dir, '/ppmi/config.R'))

#source(paste0(script_dir, '/ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))
se_rnas=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
print(prefix)

## 1. get Summarized Experiment with metrics from all time points 
## 2. Run deseq 
## 3. enrichment 

MOFAobject_clusts=MOFAobjectPD



## SETUP the script parameters here ####
# 1. select visit, 2. process mirs 
# TODO: make function to load for rnas and mirnas separately
# edit this one 
VISIT_COMP = 'V08'; VISIT = 'V08'
#VISIT=c( 'V06',  'V08') ; VISIT_COMP=VISIT

process_mirnas=FALSE ;
cell_corr=FALSE; source(paste0(script_dir,'ppmi/config.R')); outdir_s_corr0<-outdir_s

cell_corr=TRUE; source(paste0(script_dir,'ppmi/config.R')); outdir_s_corr1<-outdir_s

#cell_corr=FALSE

dir.create(outdir_s);outdir_s

if (process_mirnas){
  se_sel = se_mirs
  prefix='mirnas_'; view='miRNA'
}else{
  se_sel = se_rnas
  prefix='rnas_'; view='RNA'
}

se_filt<-filter_se(se_sel, VISIT=VISIT_COMP, sel_coh = sel_coh, sel_sub_coh = sel_ps)
table(se_filt$EVENT_ID);table(se_filt$INEXPAGE);table(se_filt$COHORT);

### Decide on the parameters settings 
## Deseq here
# Set the outdirectory 

# 1. DE files 
# 2. Venns 
# 3. Volcano plot
# 4. Enrichment analysis 
de_file
de_file = paste0(outdir_s, '/results_df.csv')
#min.count = 10
(file.exists(de_file))
if (file.exists(de_file)){
  
  
  # if de file exists load it - unfiltered de results file
  deseq2ResDF<-read.csv(paste0(de_file), row.names=1 )
}else{
  
  deseq2ResDF = deseq_by_group(se_filt, formula_deseq, min.count=min.count)
  
  if (!process_mirnas){
    # get symbols for RNA only 
    deseq2ResDF$GENE_SYMBOL<-get_symbols_vector(gsub('\\..*', '',rownames(deseq2ResDF))) 
  }
  write.csv(deseq2ResDF, de_file, row.names=TRUE)
}

deseq_all_combined<-deseq2ResDF[deseq2ResDF$mofa_sign %in% 'Significant',] # holds the significant only
deseq_all_names_combined<- gsub('\\..*', '',rownames(deseq_all_combined)) # remove the dots from ensembl names 
deseq_all_names_combined

des
length(deseq_all_names_combined)
de_genes_ppmi<-read.csv('ppmi/de_genes', header=FALSE); 
de_genes_ppmi_r<-gsub('\\..*', '',de_genes_ppmi$V1); length(de_genes_ppmi_r)

intersect(deseq_all_names_combined,de_genes_ppmi_r)




## Volcano plot #### 
volcano_title<-paste0( as.character(cell_corr), des)
pvol<-plotVolcano(  deseq2ResDF, se_filt, title=volcano_title, xlim=c(-1.1,1.1),
                      lab=deseq2ResDF$GENE_SYMBOL)
fname<-paste0(outdir_s, '/EnhancedVolcano_edited_','.jpeg')
# fname<-paste0(deseq_params, '/Volcano_', prefix, VISIT_COMP,'_cluster_', cluster_id, '.jpeg')
  
pvol
ggsave(fname,pvol, width=9,height=12, dpi=300)


#### 2. Venn from significant in top of factor

## 1. get top of factor / or all higly variable genes input into MOFA
## 2. intersect with DE 

# intersect with the top factors 
# LOAD DESEQ? 

  deseq_params = outdir_s

deseq2ResDF$log2FoldChange
order_by_metric='log2pval';order_by_metric_s='log2p'
order_by_metric='log2FoldChange'; order_by_metric_s='log2FC'

ONT='BP'
pvalueCutoff_sig=0.05
enrich_params<-paste0(ONT, '_', order_by_metric_s)
dir.create(paste0(deseq_params, '/enr/'))

enrich_compare_path=paste0(deseq_params, '/enr/', prefix, enrich_params, 'comp')


gene_list1<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )
results_file_cluster=paste0(deseq_params, '/enr/', prefix, enrich_params)
gse1<-run_enrich_per_cluster(deseq2ResDF, results_file_cluster,N_DOT=20, N_EMAP=50 )




# TODO: add mirs enrichment 


### TODO: compare with and withour correction the p-values 
# 1. whi
###COMPARE BY VISIT 

# TODO: fix by visit 
# read all visits 
deseq



deseq_all_times<-sapply( c('V04', 'V06', 'V08'), function(VISIT){
  # TODO: load deseq_file with different visit every time 
  deseq2ResDF_time<-read.csv(paste0(deseq_params_all,'/', VISIT, '/' , prefix, 'de_cluster_', cluster_id , '.csv'), row.names=1) 

  # todo: etract gene lsit
  names(gene_list1)<-gsub('\\..*', '',names(gene_list1))
  
  
  return(gene_list1)
}
)

dir.create(paste0(deseq_params_all, '/enr/'))
enrich_compare_path=paste0(deseq_params_all, '/enr/', prefix, enrich_params, cluster_id, 'time')

enrich_compare_path
gse_compare<-compareCluster(geneClusters = list(T1=deseq_all_times[[1]],T2=deseq_all_times[[2]],T3=deseq_all_times[[3]] ), 
                            fun = "gseGO", 
                            OrgDb='org.Hs.eg.db', 
                            ont=ONT, 
                            keyType = 'ENSEMBL') 


plot_enrich_compare(gse_compare,paste0(enrich_compare_path,clust_pair_s), N_EMAP = 80)

















































































































