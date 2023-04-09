

library('GOfuncR')
require(DOSE)
library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)
library('org.Hs.eg.db')
install.packages('ggridges')
require(ggridges)


padj_T=1;log2fol_T=0.00

order_by_metric<-'log2pval'

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

if (run_mofa){
  results_file = paste0(outdir, '/enrichment/gsego_')
  gse=gse_mofa}
### to run mofa results



N=10
dp<-dotplot(gse, showCategory=N, split=".sign") + facet_grid(.~.sign)
ggsave(paste0(results_file, '_dot',  '.jpeg'), width=8, height=N*0.7)

#### EMAP PLOT 
options(ggrepel.max.overlaps = Inf)

N=200
x2 <- pairwise_termsim(gse )
N=25
p<-emapplot(x2,showCategory = N,
            layout = "nicely")
p_enrich <- p + theme(text=element_text(size=12))
p_enrich

ggsave(paste0(results_file, '_emap_', N,  '.jpeg'), width=8, height=8)


#### Ridge plot: NES shows what is at the bottom of the list

r_p<-ridgeplot(gse)
r_p

#### Gene-concept plot 
gse_x <- setReadable(gse, 'org.Hs.eg.db', 'ENSEMBL')
p1_net <- cnetplot(gse_x)
node_label<-"gene"
node_label<-"all"
node_label<-"category"

N=3
p2_net<- cnetplot(gse_x, node_label=node_label, 
                       cex_label_category = 1.2, showCategory=N)

p2_net
ggsave(paste0(results_file, '_geneconcept_', node_label, '_',N, '.jpeg'), width=8, height=8)


####Visualize go terms as an undirected acyclic graph 0
goplot(gse_x)


#### heatmap
p1 <- treeplot(x2)
p2 <- treeplot(x2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')



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

