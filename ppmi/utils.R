

## Utils 
## Summarized experiment 

### TODO: move to a utils / preprocessing file because it is used also for proteoomics
library(SummarizedExperiment)


getSummarizedExperimentFromAllVisits<-function(raw_counts_all, combined){
  #
  raw_counts_all<-raw_counts_all[,!duplicated(colnames(raw_counts_all), fromLast=TRUE)]
  combined$PATNO_EVENT_ID<-paste0(combined$PATNO, '_',combined$EVENT_ID)
  
  ### some samples do not exist in metadata so filter them out 
  ## 
  common_samples<-intersect(colnames(raw_counts_all),combined$PATNO_EVENT_ID)
  unique_s<-colnames(raw_counts_all)[!(colnames(raw_counts_all) %in% common_samples)]
  metadata_filt<-combined[match(common_samples, combined$PATNO_EVENT_ID),]
  raw_counts_filt<-raw_counts_all[,match(common_samples, colnames(raw_counts_all))]
  dim(metadata_filt)[1] ==dim(raw_counts_filt)[2]
  
  
  #subset sample names

  se=SummarizedExperiment(raw_counts_filt, colData = metadata_filt)
  
  metadata_filt$COHORT_DEFINITION
  
  return(se)
  
  
}

## Create the summarized experiment by selecting VISITS and cohorts 
filter_se<-function(se, VISIT, sel_coh, sel_sub_coh=FALSE){
  
  #' Takes the raw file with all counts
  #' Filters summarized experiment by selecting VISITS and cohorts 
  #' @param VISIT
  #' @param sel_coh
  
  ##### 2.   start filtering the experiment  to normalize as appropriate 
  ## Option 1: normalize cohort and EVENT separately!! 
  # ALSO MAKE SURE THAT they are in cohort in the conversion cohort too!!
  
  
  if (length(sel_subcoh)==1 && sel_subcoh==FALSE){
        se_filt<-se[,((se$EVENT_ID %in% VISIT) & (se$COHORT %in% sel_coh ) & (se$CONCOHORT %in% sel_coh ))]
    
  }else{
      se_visit<-se[,se$EVENT_ID %in% VISIT]

      if (1 %in% sel_coh){ 
        ids<-c(se_visit$INEXPAGE %in% sel_subcoh ) ## filter the ids in parkinsons 
      }
      if (2 %in% sel_coh){
        ids<- c(ids | se_visit$INEXPAGE %in% 'INEXHC') ## also extract controls 
        
      }
      se_filt<-se_visit[, ids]
      
      
    
  }
     
  dim(se_filt)     

  
  Sample<-colnames(se_filt)
  sample_info<-DataFrame(Sample=Sample)
  
 
  
  ##### Define
  
  ### TODO: Question: Should I input everything into the matrix to normalize? 
  ### And then filter 
  
  ### batch effect and normalization 
  # Create a separate matrix with counts only
  # Include batch information if there is any
  #sample_info$Batch <- as.factor(sample_info$Batch)
  
  
  ### DEFINE THE DESEQ OBJECT with the groups appropriately 
  se_filt$EVENT_ID=as.factor(se_filt$EVENT_ID)
  se_filt$COHORT=as.factor(se_filt$COHORT)
  se_filt$PATNO=as.factor(se_filt$PATNO)
  
  
  return(se_filt)
  
}




### METADTA 

#install.packages('eeptools' )
#library('eeptools')

#as.Date(as.character(new$STATUS_DATE), format = "MM/YY",)

# TODO: fix 
get_age_at_visit<-function(new){
  AGE_AT_VISIT<-as.numeric(gsub('.*/','',new$STATUS_DATE)) - as.numeric(gsub('.*/','',new$BIRTHDT))
  return(AGE_AT_VISIT)
  }
#x_age <- age_calc( as.Date(new$BIRTHDT),          # Convert birth to age
##                   as.Date(new$STATUS_DATE),
#                  units = "years")





#### data specifc 
get_symbols_vector<-function(ens ){
  #' @param ens ensemble ids to conver to symbols 
  #' @returns symbols_ordered the total 
  #'  
  #'  
  
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  symbols <- symbols[!is.na(symbols)]
  symbols_ordered <- symbols[match(ens, names(symbols))]
  na_ind<-is.na(symbols_ordered);
  
  # Add ensembl ids if no symbol found
  symbols_ordered[na_ind]=ens[na_ind]
  return(symbols_ordered)
  
  
}



######## DE ANALYSIS #######
#results_de<-mark_signficant(
  # test
  #de_res= results_de
  #padj_T = padj_T_overall; log2fol_T = log2fol_T_overall; padj_name ='adj.P.Val'
  #log2fc_name = 'logFC'   
  #outdir_single=outdir_s_p

mark_signficant<-function(de_res, padj_T, log2fol_T, padj_name='padj', log2fc_name='log2FoldChange', outdir_single=outdir_s ){
  ## mark a significant column and write to file
  
  signif_file<-paste0(outdir_single,'/significant', padj_T, '_',log2fol_T, '.csv')
  
  de_res$significant <- ifelse(de_res[,padj_name] < padj_T , "Significant", NA)
  de_res$sign_lfc <- ifelse(de_res[,padj_name] < padj_T & abs(de_res[,log2fc_name]) >log2fol_T , "Significant", NA)
  # Examine this data frame
  # Order the significant to save as a new output file 
  head(de_res)
  # LARGER ONE not saved 
  sign_only<-de_res[de_res$sign_lfc=='Significant',]
  sign_only_ordered<-sign_only[order(sign_only[,padj_name], decreasing = FALSE),]
  sign_only_ordered<-sign_only[order(-sign_only[,'abslog2pval'], decreasing = FALSE),]
  
  sign_only_ordered<-sign_only_ordered[!is.na(de_res$sign_lfc),]
  write.csv(sign_only_ordered,signif_file, row.names = TRUE)
  ### create also a more strict file? 
  return(de_res)
}



######## ENRICHMENT ANALYSIS 


get_ordered_gene_list<-function(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 ){
  
  #### Gives a gene list cut by the thresholds padj_T and log2fol and orders by specific metric supplied 
  #' @param padj_T filter the genes by metric padj_T
  #' @param log2fol_T  filter the genes by metric  log2fol_T, default: 0
  #' @param order_by_metric metric to order the gene list by 
  res=deseq2ResDF
  #res=deseq2ResDF_2
  
  res$sign_lfc <- ifelse(res$padj <padj_T & abs(res$log2FoldChange) >log2fol_T , "Significant", NA)
  
  length(which(!is.na(res$sign_lfc )))
  res=res[res$sign_lfc=='Significant'& !is.na(res$sign_lfc),]
  
  
  # Order the DE gene list by the stat statistic 
  #remove negatives thatw ere introduced with vst transofrmations
  
  res$log2pval<-res$log2FoldChange*-log10(res$padj)
  res$signlog2pval<-sign(res$log2FoldChange)*-log10(res$padj)
  res<-res[res$baseMean>0,]
  
  #res <- res[order(-res$stat),]
  res <- res[order(-res[,order_by_metric]),]
  gene_list<-res[, order_by_metric]
  names(gene_list)<-rownames(res)
  
  return(gene_list)
  
  
  
}

write_filter_gse_results<-function(gse_full,results_file,pvalueCutoff, pvalueCutoff_sig=0.05  ){
  
  ### Takes the full gse results, ie. without threshold significance, 
  # saves it, 
  # filters it by pvalueCutoff_sig
  # and saves the filter 
  #' @param gse_full full gse results objects 
  #' @param results_file the file name to write results  (without .csv)
  #' @param pvalueCutoff the pvalue used to obtain the gse results 
  pval_to_use='p.adjust'
  write.csv(as.data.frame(gse_full@result), paste0(results_file, pvalueCutoff, '.csv'))
  gse_sig_result<-gse_full@result[gse_full@result[,pval_to_use]<pvalueCutoff_sig,]
  write.csv(as.data.frame(gse_sig_result), paste0(results_file, pvalueCutoff_sig, '.csv'))
  
  # rewrite
  dim(gse_full); dim(gse_sig_result)
  ## filter gse result to significant only 
  gse=dplyr::filter(gse_full, p.adjust < pvalueCutoff_sig)
  return(gse)
}


run_enrichment_plots<-function(gse, results_file,N_EMAP=25, N_DOT=15, N_TREE=30, N_NET=30, showCategory_list=FALSE, process_mofa=FALSE){
  
  require(clusterProfiler)
  
  require('GOfuncR')
  require(DOSE)
  
  
  
  N=25
  ## TODO: ADD FACET IF SIGNED 
  if (length(showCategory_list)>1){
    print('Filter by selected category')
    N_DOT<-showCategory_list
    N_TREE<-showCategory_list
    N_NET<-showCategory_list
    N_EMAP<-showCategory_list
    
    write_n=FALSE
  }
  
  ### print a signed and unsigned dotplot 
  # because it does not make sense if we dont rank by logFC
  # or in the mofa case where we rank by importance in factor 
  
  dp<-dotplot(gse, showCategory=N_DOT, 
              font.size=15
  )
  dp<-dp+theme(axis.ticks=element_blank() , 
               axis.text.x = element_blank())
  show(dp)
  
  
  
  if (process_mirnas){
    width=6}else{width=6}
  
  ggsave(paste0(results_file, '_dot', N_DOT, '.jpeg'), 
         plot=dp, width=width, height=N_DOT*0.5, 
         dpi = 300)
  
  if (!(process_mirnas) & !(run_ORA)){
    
    dp_sign<-dotplot(gse, showCategory=N_DOT, split=".sign") + facet_grid(.~.sign)
    ggsave(paste0(results_file, '_dot_sign', N_DOT,  '.jpeg'), width=8, height=N*0.7)
    
  }
  
  #### EMAP PLOT 
  #N_EMAP=50
  options(ggrepel.max.overlaps = Inf)
  x2 <- pairwise_termsim(gse )
  #if (process_mirnas){N=15}
  p<-emapplot(x2,showCategory = N_EMAP,
              layout = "nicely", 
              cex_label_category=0.8)
  p_enrich <- p + theme(text=element_text(size=12))
  p_enrich
  
  if (is.numeric(N_EMAP)){write_n=N_EMAP}
  ggsave(paste0(results_file, '_emap_', write_n,  '.jpeg'), width=9, height=9, 
         dpi = 300)
  
  
  #### Ridge plot: NES shows what is at the bottom of the list
  
  
  N_RIDGE=25
  # only if all 3 are false run it 
  if ( !(process_mirnas) && !(process_mofa) && !(run_ORA)){
    print('ridge')
    r_p<-ridgeplot(gse, showCategory = N_RIDGE)
    r_p
    ggsave(paste0(results_file, '_ridge_', N_RIDGE, '.jpeg'), width=8, height=8)
  }
  
  
  
  
  #### Gene-concept plot 
  
  
  
  if (gse@keytype=='ENSEMBL' ){
    gse_x <- setReadable(gse, 'org.Hs.eg.db', 'ENSEMBL')
    
  }else{
    gse_x=gse
    
  }
  
  
  p1_net <- cnetplot(gse_x)
  
  node_label<-"gene"
  node_label<-"category"
  node_label<-"all"
  
  p2_net<- cnetplot(gse_x,
                    node_label=node_label,
                    cex_label_category = 1.2, showCategory=N_NET)
  
  p2_net
  if (is.numeric(N_NET)){write_n=N_NET}
  
  ggsave(paste0(results_file, '_geneconcept_', node_label, '_',write_n, '.jpeg'), width=8, height=8)
  
  
  ####Visualize go terms as an undirected acyclic graph 0
  if (!process_mirnas){
    goplot(gse_x)
    ggsave(paste0(results_file, '_goplot_', node_label, '_',write_n, '.jpeg'), width=8, height=8)
    
  } 
  library(ggtree)
  library(ggplot2)
  
  #install.packages('ggtree')
  
  #### heatmap
  N_TREE=16
  p1 <- treeplot(x2,showCategory =N_TREE, nWords=0)
  p1
  p2_tree <- treeplot(x2, hclust_method = "average", 
                      showCategory =N_TREE, nWords=0, 
                      #offset_tiplab=5, 
                      label_format =50, 
                      fontsize = 300, 
                      extend=-0.001, 
                      offset=15, 
                      hilight=FALSE, 
                      branch.length=0.1)
  
  #aplot::plot_list(p1, p2_tree, tag_levels='A')
  #ggsave(paste0(results_file, '_clusterplot_', node_label, '_',N, '.jpeg'), width=8, height=8)
  
  p2_tree
  #write_n='test'
  ggsave(paste0(results_file, '_clusterplot_average_',write_n, '.jpeg'),
         width=10, height=0.4*N_TREE, dpi=300)
  
  
  return(list(dp, p_enrich, p2_tree))
  
}





#### Configuration 

VISIT='V08'

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
VISIT='V08'
source(paste0(script_dir, '/config.R'))
gene_list<-get_genelist_byVisit(VISIT)


source(paste0(script_dir, '/config.R'))
outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)

### setup
length(gene_list)
tail(gene_list)
ONT='BP'



results_file_rnas<-paste0(outdir_enrich, '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)
res_path_rnas<-paste0(results_file_rnas, 'gse.RDS')

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


