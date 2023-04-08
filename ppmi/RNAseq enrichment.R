
BiocManager::install('clusterProfiler', force=TRUE)
#BiocManager::install('apeglm')
#install.packages('ggnewscale')
#BiocManager::install('GOfuncR')
#install.packages("BiocManager")

library('GOfuncR')
require(DOSE)
library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)

padj_T=1;log2fol_T=0.00

order_by_metric<-'log2pval'

deseq2Results = read.csv(paste0(outdir_s, '/results.csv'), row.names = 1)
deseq2ResDF <- as.data.frame(deseq2Results) 

 
gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )


names(gene_list)<-gsub('\\..*', '',names(gene_list))
length(gene_list)

ONT='BP'

outdir_enrich<-paste0(outdir_s,'/enrichment/')
dir.create(outdir_enrich)
                      
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



require(DOSE)
library('enrichplot')
#BiocManager::install('ggnewscale')
#library('ggnewscale')
gse@result

write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))
N=10
dp<-dotplot(gse, showCategory=N, split=".sign") + facet_grid(.~.sign)
ggsave(paste0(results_file, '_dot',  '.jpeg'), width=8, height=N*0.7)

#### EMAP PLOT 
options(ggrepel.max.overlaps = Inf)


x2 <- pairwise_termsim(gse)
N=25
p<-emapplot(x2,showCategory = N,
            layout = "nicely")
p_enrich <- p + theme(text=element_text(size=12))
p_enrich

ggsave(paste0(results_file, '_emap_', N,  '.jpeg'), width=8, height=8)


#### Gene-concept plot 
gse_x <- setReadable(gse, 'org.Hs.eg.db', 'ENSEMBL')
p1 <- cnetplot(gse_x, gene_list)
p1

####Visualize go terms as an undirected acyclic graph 0
goplot(gse_x)



############# KEGG


parents<-get_parent_nodes(gse_x$ID, term_df = NULL, graph_path_df = NULL, godir = NULL)

parents$distance


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

