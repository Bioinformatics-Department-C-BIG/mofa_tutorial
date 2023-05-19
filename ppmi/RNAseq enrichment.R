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





## Filter 
gse=write_filter_gse_results(gse_full, results_file_rnas, pvalueCutoff)

#### 
# run all results

require(DOSE)
library('enrichplot')


#results_file=mir_results_file_anticor
#gse = gse_mirnas
pvalueCutoff_sig=0.01
sel<-gse_full@result$pvalue<pvalueCutoff_sig
gse=filter(gse_full, p.adjust < pvalueCutoff_sig)


enrich_plots<-run_enrichment_plots(gse=gse,results_file=results_file_rnas, N_DOT=15, N_EMAP=25 )
dp=enrich_plots[[1]]
p_enrich=enrich_plots[[2]]
p2_tree=enrich_plots[[3]]



#### 
gse@result[gse@result$pvalue>0.05,]

df<-gse@result
hist_p<-ggplot(df, aes(x=-log10(pvalue))) + 
  geom_histogram(color="black", fill="white")
hist_p

ggsave(paste0(results_file_rnas, '_pval_hist.png'),plot=last_plot() )




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

