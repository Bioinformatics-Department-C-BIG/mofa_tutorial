#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
#source(paste0(script_dir,'ppmi/setup_os.R'))


#detach("package:AnnotationDbi", unload=TRUE)
#detach("package:org.Hs.eg.db", unload=TRUE)

library('GOfuncR')
require(DOSE)
#library('apeglm')
library(clusterProfiler)
#BiocManager::install('clusterProfiler')
#library(AnnotationDbi)
#library(ensembldb)
#install.packages('ggridges')
require(ggridges)
library(ggplot2)

suppressWarnings(library('R.filesets' ))
library('enrichplot' )





#### Configuration 

#VISIT='V08'
# Rerun config 
# Load Visit from the config and not from here? 
source(paste0(script_dir, 'ppmi/config.R'))


process_mirnas<-FALSE
padj_T=1;log2fol_T=0.00;order_by_metric<-'log2pval'


get_genelist_byVisit<-function(VISIT){
  
  ### Input visit AND return list 
  ##'
  ##'
  
  deseq2ResDF_2 = as.data.frame(read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1))
  gene_list<-get_ordered_gene_list(deseq2ResDF_2,  order_by_metric, padj_T=1, log2fol_T=0 )
  names(gene_list)<-gsub('\\..*', '',names(gene_list))
  return(gene_list)
}
#source(paste0(script_dir, '/config.R'))
gene_list<-get_genelist_byVisit(VISIT)

outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)

### setup
length(gene_list)
tail(gene_list)
ONT='BP'



results_file<-paste0(outdir_enrich, '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)
res_path<-paste0(results_file, 'gse.RDS')

#### Run and return the whole set of p-values with pcutoff=1 
## then filter 
if (file.exists(res_path)){
  gse_full=loadRDS(res_path)
  
}else{
  
  pvalueCutoff<-1
  gse_full <- clusterProfiler::gseGO(gene_list, 
                                     ont=ONT, 
                                     keyType = 'ENSEMBL', 
                                     OrgDb = 'org.Hs.eg.db', 
                                     pvalueCutoff  = pvalueCutoff)
  saveRDS(gse_full, res_path)
  
}

## Filter 
gse=write_filter_gse_results(gse_full, results_file, pvalueCutoff)

#### 
# run all results

require(DOSE)
library('enrichplot')

results_file=results_file
gse = gse

#results_file=mir_results_file_anticor
#gse = gse_mirnas
pvalueCutoff_sig=0.05
sel<-gse_full@result$pvalue<pvalueCutoff_sig
gse=filter(gse_full, p.adjust < pvalueCutoff_sig)

#text_p<-get_pval_text(gse_full, pvalueCutoff_sig)
text_p<-''

enrich_plots<-run_enrichment_plots(gse=gse,results_file=results_file, N_DOT=15, N_EMAP=25, text_p=text_p )
dp=enrich_plots[[1]]
p_enrich=enrich_plots[[2]]
p2_tree=enrich_plots[[3]]




