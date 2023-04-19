
#install.packages('VennDiagram')
library(rbioapi)
library('VennDiagram')
library(enrichplot)
VISIT='BL'
process_mirnas<-TRUE
source(paste0(script_dir, '/config.R'))



order_by_metric<-'abslog2pval'
order_by_metric<-'abslog2pval'
order_by_metric<-'log2pval'
order_by_metric<-'abslog2pval'
order_by_metric<-'log2FoldChange'
order_by_metric<-'log2pval'
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

Padj_T_paths=0.01

deseq2ResDF = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)
outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)

gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T, log2fol_T )
top_n=length(gene_list);
#top_n=200
top_n
gene_list_cut<-gene_list[1:top_n]
length(gene_list); length(gene_list_cut)
mirs=names(gene_list_cut)
mirs
enrich_params<-paste0('_', padj_T, '_',  log2fol_T, '_',  order_by_metric)
mir_results_file<-paste0(outdir_enrich, '/mirs_enrich_', enrich_params)

#### Up to here we can take the output for other methods 




############## RUN MIEAA ######################################

gsea_results_fname<-paste0(mir_results_file,'_mieaa_res.csv' )


if (file.exists(gsea_results_fname)){
  ### Load enrichment results if available
  mieaa_all_gsea<-read.csv(gsea_results_fname, header=TRUE)
  
}else{
  ## otherwise run GSEA analysis 
  mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                                     mirna_type = "mature",
                                     test_type = "GSEA",
                                     species = 'Homo sapiens',
                                     sig_level=pvalueCutoff
  )
  ### write to file to load next time 
  write.csv(mieaa_all_gsea, gsea_results_fname, row.names = FALSE)
  
}


## remove . and \ from mir names 
colnames(mieaa_all_gsea)<-gsub('-','.', colnames(mieaa_all_gsea))
colnames(mieaa_all_gsea)<-gsub('/','.', colnames(mieaa_all_gsea))




#table(mieaa_all_gsea$Category)
Category<-'Annotation (Gene Ontology)';
Category<-'Gene Ontology (miRWalk)';
Category<-'GO Biological process (miRPathDB)';


mieaa_gsea_1<-mieaa_all_gsea[mieaa_all_gsea$Category==Category,]

write.csv(mieaa_gsea_1, paste0(mir_results_file, '_', pvalueCutoff,'.csv' ))
hist(mieaa_gsea_1$P.adjusted)

mieaa_gsea_1_cut<-mieaa_gsea_1[mieaa_gsea_1$P.adjusted<Padj_T_paths, ]
mir_paths<-mieaa_gsea_1_cut[,c(2)]


results_df<-paste0(mir_results_file, '_', Category)

write.csv(mieaa_all_gsea, paste0(mir_results_file, '_', '.csv' ))




##### plot results 
##
#DOSE::enrichResult
#  EnrichResult(mieaa_all_gsea)

  

df=mieaa_gsea_1
df$P.adjusted<-as.numeric(df$P.adjusted)
df$padj<-as.numeric(df$P.adjusted)


df$Observed<-as.numeric(df$Observed)

df_ord<-df[order(df$padj),]

df_ord<-df_ord[1:40,]
#df_ord$padj

mir_enrich_p<-ggplot(df_ord, aes(x=reorder(Subcategory, padj), y=Observed, fill=padj))+
  geom_bar(position='dodge', stat='identity')+
  coord_flip()



ggsave(paste0(mir_results_file, '_', Category, '_bar',  '.png'), height = 7, width=8)

dir.create(paste0(outdir_enrich, Category))

mir_results_file_by_cat<-paste0(outdir_enrich, Category, '/mirs_enrich_', enrich_params)


######### convert to enrichResult to use gsego functios


#install.packages("remotes")
#remotes::install_github("jmw86069/jamenrich")
library('multienrichjam')
library('clusterProfiler')

mieaa_gsea_1$P.adjusted<-as.numeric(mieaa_gsea_1$P.adjusted)


mieaa_gsea_1$keyColname=mieaa_gsea_1$Subcategory
mieaa_gsea_1_ord=mieaa_gsea_1[order(mieaa_gsea_1$P.adjusted),]
#mieaa_gsea_1_ord_prob<-mieaa_gsea_1_ord
enr_full <- multienrichjam::enrichDF2enrichResult(as.data.frame(mieaa_gsea_1_ord),
                                             keyColname =  'Subcategory',
                                             geneColname ='miRNAs.precursors',
                                             pvalueColname = 'P.adjusted', 
                                             pvalueCutoff = pvalueCutoff)

# descriptionColname = "Subcategory",

## double check why they are CUT 

results_file=mir_results_file_by_cat
mir_results_file_by_cat
gse_sig=write_filter_gse_results(gse_full=enr_full, mir_results_file_by_cat, pvalueCutoff)
paste0(mir_results_file_by_cat, pvalueCutoff, '.csv')

mir_results_file_by_cat

#write.csv(gse_sig@result, paste0(mir_results_file_by_cat,pvalueCutoff, '.csv'))



### requires source('RNAseq enrichment.R') # TODO: MOVE TO A UTILS SCRIPT 
run_enrichment_plots(gse=enr, results_file=mir_results_file_by_cat, N_EMAP=15)


############
###########



