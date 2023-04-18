

library('GOfuncR')
require(DOSE)
library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)
library('org.Hs.eg.db')
#install.packages('ggridges')
require(ggridges)
suppressWarnings(library('R.filesets' ))
library('enrichplot' )

VISIT='V08'
process_mirnas<-FALSE
source(paste0(script_dir, '/config.R'))



padj_T=1;log2fol_T=0.00;order_by_metric<-'log2pval'

deseq2Results = read.csv(paste0(outdir_s, '/results.csv'), row.names = 1)
deseq2ResDF <- as.data.frame(deseq2Results) 
outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)

### setup


gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )
names(gene_list)<-gsub('\\..*', '',names(gene_list))
length(gene_list)

ONT='BP'



results_file<-paste0(outdir_enrich, '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)

res_path<-paste0(results_file, 'gse.RDS')

if (file.exists(res_path)){
  gse=loadRDS(res_path)
  
}else{
  
  
  gse <- clusterProfiler::gseGO(gene_list, 
                                ont=ONT, 
                                keyType = 'ENSEMBL', 
                                OrgDb = 'org.Hs.eg.db', 
                                pvalueCutoff  = 0.05)
  saveRDS(gse, res_path)
  
}

write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))



#### 
# run all results

require(DOSE)
library('enrichplot')

results_file=results_file
gse = gse

results_file=mir_results_file_anticor
gse = gse_mirnas



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
  if (process_mirnas){
    
    #results_file=mir_results_file_anticor
    #gse=gse_mirnas;
    
    dp<-dotplot(gse, showCategory=N_DOT)
    
  }else{
    dp<-dotplot(gse, showCategory=N_DOT, split=".sign") + facet_grid(.~.sign)
    
  }
  
  ggsave(paste0(results_file, '_dot',  '.jpeg'), width=8, height=N*0.7)
  
  #### EMAP PLOT 
  options(ggrepel.max.overlaps = Inf)
  
  x2 <- pairwise_termsim(gse )
  x2
  #if (process_mirnas){N=15}
  p<-emapplot(x2,showCategory = N_EMAP,
              layout = "nicely")
  p_enrich <- p + theme(text=element_text(size=12))
  p_enrich
  
  if (is.numeric(N_EMAP)){write_n=N_EMAP}
  ggsave(paste0(results_file, '_emap_', write_n,  '.jpeg'), width=9, height=9)
  
  
  #### Ridge plot: NES shows what is at the bottom of the list
  
  if (!process_mirnas){
    r_p<-ridgeplot(gse)
    r_p
    ggsave(paste0(results_file, '_ridge_.jpeg'), width=8, height=8)
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
  
  
  #### heatmap
  p1 <- treeplot(x2,showCategory =N_TREE)
  p2_tree <- treeplot(x2, hclust_method = "average", showCategory =N_TREE )
  #aplot::plot_list(p1, p2_tree, tag_levels='A')
  #ggsave(paste0(results_file, '_clusterplot_', node_label, '_',N, '.jpeg'), width=8, height=8)
  
  p2_tree
  ggsave(paste0(results_file, '_clusterplot_average_',write_n, '.jpeg'), width=12, height=8)
  
  
  return(list(dp, p_enrich, p2_tree))
  
}



enrich_plots<-run_enrichment_plots(gse=gse,results_file=results_file )
dp=enrich_plots[[1]]
p_enrich=enrich_plots[[2]]
p2_tree=enrich_plots[[3]]


#### RUN ENRICHMENT WITH FILTERED PATHS FROM 3 MODALITIES ####

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

