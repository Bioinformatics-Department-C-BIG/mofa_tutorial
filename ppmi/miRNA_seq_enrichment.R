#source('ppmi/setup_os.R')


#install.packages('VennDiagram')
library(rbioapi)
library('VennDiagram')
library(enrichplot)
process_mirnas<-TRUE
source(paste0(script_dir, 'ppmi/config.R'))
source(paste0(script_dir, 'ppmi/utils.R'))
library(ggplot2)





order_by_metric<-'abslog2pval'
order_by_metric<-'abslog2pval'
order_by_metric<-'log2pval'
order_by_metric<-'abslog2pval'
order_by_metric<-'log2pval'
order_by_metric<-'log2FoldChange'

VISIT
use_anticor=FALSE
## UPDATED TO 0,1--DO not filter the list of de mirnas, since it runs gsea it needs the complete list 
if (VISIT=='V08'){

  padj_T=0.01
  log2fol_T=0.0
  padj_T=0.05
  log2fol_T=0.1
  ## MORE STRICT IN V08
  padj_T=0.01
  log2fol_T=0.1
  padj_T=1
  log2fol_T=0
  
}else{
  padj_T=1
  log2fol_T=0
}
  
if (use_anticor){
  order_by_metric<-'log2pval'
  padj_T=0.05
  log2fol_T=0
  
}

log2fol_T;padj_T;
log2fol_T
### run the enrichment if it has not been ran yet!! 
padj_T;log2fol_T
Padj_T_paths=0.01

deseq2ResDF = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)
de_mirs<-deseq2ResDF
outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)



gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T, log2fol_T )
gene_list_sig_only<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=0.05, log2fol_T=0.1 )

test_type='GSEA';top_n=length(gene_list);
test_type='ORA'; top_n=100
### BL: there is nothing significant to plot... whyyyy

gene_list_cut_ora<-gene_list_sig_only[order(abs(gene_list_sig_only))]
length(gene_list_cut_ora)
gene_list_cut_ora

gene_list_cut<-gene_list[order(abs(gene_list))][1:top_n]
gene_list_cut<-gene_list_cut[order(gene_list_cut, decreasing=TRUE)]

if (test_type=='ORA'){
  gene_list_cut=gene_list_cut_ora
  mirs_ora=names(gene_list_cut_ora)
}

length(gene_list); length(gene_list_cut)
mirs=names(gene_list_cut)
mirs

enrich_params<-paste0('_', padj_T, '_',  log2fol_T, '_',  order_by_metric,test_type,  top_n)
mir_results_file<-paste0(outdir_enrich, '/mirs_enrich_', enrich_params)

#### Up to here we can take the output for other methods 

gene_list_cut


############## RUN MIEAA ######################################

gsea_results_fname<-paste0(mir_results_file,'_mieaa_res.csv' )
pvalueCutoff=1
mirs=gsub( '\\.','-', selected_mirs)
mirs=gsub( '\\.','-', de_group_vs_control_and_time2)
mirs
test_type='ORA'
if (file.exists(gsea_results_fname)){
  ### Load enrichment results if available
  mieaa_all_gsea<-read.csv(gsea_results_fname, header=TRUE)
  ### TODO: Rerun with updated pvalue cutoff 
}else{
  ## otherwise run GSEA analysis 
  
  
  mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                                     mirna_type = "mature",
                                     test_type = test_type,
                                     species = 'Homo sapiens',
                                     sig_level=pvalueCutoff
  )


  write.csv(mieaa_all_gsea, gsea_results_fname, row.names = FALSE)
  
}

colnames(mieaa_all_gsea)<-make.names(  colnames(mieaa_all_gsea))
mieaa_all_gsea_sig<-mieaa_all_gsea %>%
  dplyr::filter(Category %in%c('GO Biological process (miRPathDB)')) %>%
    dplyr::filter(P.adjusted<0.05)

mieaa_all_gsea_sig$Subcategory
View(mieaa_all_gsea_sig)

mieaa_targets<-mieaa_all_gsea %>%
  dplyr::filter(Category %in%c('Target genes (miRTarBase)')) %>%
  dplyr::filter(P.adjusted<0.05)

View(mieaa_targets)
### write to file to load next time 

### TODO: CREATE FUNCTION OF ENRICHMENT 
### FROM HERE ONWARDS RUN MOFA 
## TODO: MAKE THIS A FUNCTION that takes in arguments of mieaa_all_gsea and mir_results_file 
process_mofa=FALSE
if (process_mofa){
  
  mieaa_all_gsea=mieaa_all_gsea_mofa
  mir_results_file=paste0(outdir, '/enrichment/mirs_')
  outdir_enrich=paste0(outdir, '/enrichment/', factor)
}


## remove . and \ from mir names 


use_mofa=FALSE
if (use_mofa){
  for (fn in 1:nfactors){
    gse_mofa_mirs=list_mirs[[fn]]
    
    mieaa_all_gsea=gse_mofa_mirs
    
    mieaa_res<-mirna_enrich_res_postprocessing(mieaa_all_gsea, Category)
    mieaa_gsea_1=mieaa_res[[1]]
    enr_full=mieaa_res[[2]]
    print(any(enr_full@result$p.adjust<0.05))
  }
}






dir.create(paste0(outdir_enrich, Category))


mir_results_file_by_cat<-paste0(outdir_enrich, Category, '/mirs_enrich_', enrich_params)
######### convert to enrichResult to use gsego functios
results_file=mir_results_file_by_cat
mir_results_file_by_cat
### Posto process and return enrichresult 
mieaa_res<-mirna_enrich_res_postprocessing(mieaa_all_gsea, Category = Category, mir_results_file = mir_results_file_by_cat)
mieaa_gsea_1=mieaa_res[[1]]
enr_full=mieaa_res[[2]]
enr_full

##### plot results 
##
#DOSE::enrichResult
#  EnrichResult(mieaa_all_gsea)
### intermediate plot --- remove..??


df=mieaa_gsea_1
min(mieaa_gsea_1$P.adjusted) ### with LOGFC IT WORKS BETTER

df$P.adjusted<-as.numeric(df$P.adjusted)
df$padj<-as.numeric(df$P.adjusted)
df$Observed<-as.numeric(df$Observed)
df_ord<-df[order(df$padj),]
df_ord<-df_ord[1:40,]

mir_enrich_p<-ggplot(df_ord, aes(x=reorder(Subcategory, padj), y=Observed, fill=padj))+
  geom_bar(position='dodge', stat='identity')+
  coord_flip()
ggsave(paste0(mir_results_file, '_', Category, '_bar',  '.png'), plot=mir_enrich_p, height = 7, width=8)




gse_sig=write_filter_gse_results(gse_full=enr_full, mir_results_file_by_cat, pvalueCutoff_sig)

gse_sig@result$p.adjust

paste0(mir_results_file_by_cat, pvalueCutoff, '.csv')

mir_results_file_by_cat

#write.csv(gse_sig@result, paste0(mir_results_file_by_cat,pvalueCutoff, '.csv'))

#text_p<-get_pval_text(enr_full, pvalueCutoff_sig)
text_p=''
process_mirnas
### requires source('RNAseq enrichment.R') # TODO: MOVE TO A UTILS SCRIPT 
run_enrichment_plots(gse=gse_sig, results_file=mir_results_file_by_cat, N_EMAP=15, N_DOT=15, text_p=text_p)
dim(gse_sig)

############
###########



