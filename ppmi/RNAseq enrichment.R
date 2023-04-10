




  
require(DOSE)


library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)

padj_T=1;log2fol_T=0.00

order_by_metric<-'log2pval'


 
gene_list<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )


names(gene_list)<-gsub('\\..*', '',names(gene_list))
length(gene_list)

ONT='BP'

gse <- clusterProfiler::gseGO(gene_list, 
                              ont=ONT, 
                              keyType = 'ENSEMBL', 
                              OrgDb = 'org.Hs.eg.db', 
                              pvalueCutoff  = 0.05)





require(DOSE)
library('enrichplot')
#BiocManager::install('ggnewscale')
library('ggnewscale')
gse@result

results_file<-paste0(outdir_s, '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T, order_by_metric)
write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))
N=10
dp<-dotplot(gse, showCategory=N, split=".sign") + facet_grid(.~.sign)
ggsave(paste0(results_file, '_dot',  '.jpeg'), width=8, height=N*0.7)

#### EMAP PLOT 

x2 <- pairwise_termsim(gse)
N=25
p<-emapplot(x2,showCategory = N,
            layout = "nicely")
p_enrich <- p + theme(text=element_text(size=12))


ggsave(paste0(results_file, '_emap_', N,  '.jpeg'), width=8, height=8)




############# KEGG






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

