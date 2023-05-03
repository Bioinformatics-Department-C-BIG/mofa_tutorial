script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir,'/setup_os.R'))


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




write_filter_gse_results<-function(gse_full,results_file,pvalueCutoff  ){
  
  ### Takes the full gse results, ie. without threshold significance, 
  # saves it, 
  # filters it by pvalueCutoff_sig
  # and saves the filter 
  #' @param gse_full full gse results objects 
  #' @param results_file the file name to write results  (without .csv)
  #' @param pvalueCutoff the pvalue used to obtain the gse results 
  pval_to_use='p.adjust'
  write.csv(as.data.frame(gse_full@result), paste0(results_file, pvalueCutoff, '.csv'))
  pvalueCutoff_sig<-0.05
  gse_sig_result<-gse_full@result[gse_full@result[,pval_to_use]<pvalueCutoff_sig,]
  write.csv(as.data.frame(gse_sig_result), paste0(results_file, pvalueCutoff_sig, '.csv'))
  
  # rewrite
  dim(gse_full); dim(gse_sig_result)
  ## filter gse result to significant only 
  gse=dplyr::filter(gse_full, p.adjust < pvalueCutoff_sig)
  return(gse)
}


run_enrichment_plots<-function(gse, results_file,N_EMAP=25, N_DOT=15, N_TREE=30, N_NET=30, showCategory_list=FALSE){
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
    width=4}else{width=6}
  
  ggsave(paste0(results_file, '_dot', N_DOT, '.jpeg'), 
         plot=dp, width=width, height=N_DOT*0.5, 
         dpi = 300)
  
  if (!(process_mirnas || run_ORA)){
    
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
  
  if (!process_mirnas & !process_mofa){
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
pvalueCutoff_sig=0.01
sel<-gse_full@result$pvalue<pvalueCutoff_sig
gse=filter(gse_full, p.adjust < pvalueCutoff_sig)


enrich_plots<-run_enrichment_plots(gse=gse,results_file=results_file, N_DOT=15, N_EMAP=25 )
dp=enrich_plots[[1]]
p_enrich=enrich_plots[[2]]
p2_tree=enrich_plots[[3]]


results_file

#### 
gse@result[gse@result$pvalue>0.05,]

df<-gse@result
hist_p<-ggplot(df, aes(x=-log10(pvalue))) + 
  geom_histogram(color="black", fill="white")
hist_p

ggsave(paste0(results_file, '_pval_hist.png'),plot=last_plot() )




#### RUN ENRICHMENT WITH FILTERED PATHS FROM 3 MODALITIES ####
run_anova=FALSE
use_pval=TRUE
### prerequisite: source venn diagrams
out_compare<-'ppmi/plots/single/compare/'
int_params<-paste0(padj_paths, '_', VISIT, '_p_anova_',run_anova, 'pval_', use_pval )

intersection_all_three_plot<-read.csv(paste0(out_compare,'interesction_pathways' , int_params, '.csv') )
#as.character(intersection_all_three_plot)
showCategory_list=as.vector(unlist(intersection_all_three_plot))
length(showCategory_list)


gse_common<-gse %>% 
  dplyr::filter(Description %in%showCategory_list )


enrich_plots<-run_enrichment_plots(gse=gse_common,
                                   results_file=paste0(results_file, '_intersect_',int_params ),
                                   N_EMAP=20, N_NET=20)


write.csv(gse_common@result, paste0(results_file, '_intersect_',int_params, '.csv' ))
run_mofa=FALSE
#run_mofa=TRUE
### TODO: FIX AND MAKE A FUNCTION OF THIS SO I CAN USE IN MOFA TOO 
if (run_mofa){
  for (factor in c(1:8)){
    
    
    results_file_mofa = paste0(outdir, '/enrichment/gsego_',factor,'_')
    gse_mofa=list1[[factor]]
    write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))
    
    ### to run mofa results
    run_enrichment_plots(gse=gse_mofa, results_file = results_file_mofa)
    
    
    if (process_mirnas){
      results_file=mir_results_file_by_cat
      gse=enr;
      
      
    }
    
    
    
    ############# KEGG
    
    
    #parents<-get_parent_nodes(gse_x$ID, term_df = NULL, graph_path_df = NULL, godir = NULL)
    
    #parents$distance
    
  }
}


### CLUSTER COMPARE



#### 
#library(ggplot2)
#library(dplyr)
#library(stringr)
#gse@result
#x=gse
### count the gene number
#gene_count<- x@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
#
### merge with the original dataframe
#dot_df<- left_join(x@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
#
### plot
#library(forcats) ## for reordering the factor
#ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
#  geom_point(aes(size = GeneRatio, color = p.adjust)) +
#  theme_bw(base_size = 14) +
#  scale_colour_gradient(limits=c(0, 0.10), low="red") +
#  ylab(NULL) +
#  ggtitle("GO pathway enrichment")
#
#
#
#
#
#
#gsea_run <- clusterProfiler::GSEA(gene_list, 
#                              keyType = 'ALIAS', 
#                              pAdjustMethod = 'BH')
#
#
#

