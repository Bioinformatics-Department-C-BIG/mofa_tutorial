
#BiocManager::install('clusterProfiler', force=TRUE)
#BiocManager::install('apeglm')
# BiocManager::install('AnnotationDbi')


  
require(DOSE)


library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)

res=deseq2ResDF

log2fol_T=0.1
padj_T=0.05
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < padj_T , "Significant", NA)
deseq2ResDF$sign_lfc <- ifelse(deseq2ResDF$padj < padj_T & abs(deseq2ResDF$log2FoldChange) >T_lfc , "Significant", NA)




padj_T=0.05;log2fol_T=0.1
res=deseq2ResDF
res$sign_lfc <- ifelse(res$padj <padj_T & abs(res$log2FoldChange) >log2fol_T , "Significant", NA)

length(which(!is.na(res$sign_lfc )))


res=res[res$sign_lfc=='Significant'& !is.na(res$sign_lfc),]

dim(res)
res
# Order the DE gene list by the stat statistic 
#remove negatives thatw ere introduced with vst transofrmations


res<-res[res$baseMean>0,]

#res <- res[order(-res$stat),]
res <- res[order(-res$log2FoldChange),]

gene_list<-res[,'log2FoldChange']

names(gene_list)<-rownames(res)


#  Takes input the DE genes from DESeq2



names(gene_list)<-gsub('\\..*', '',names(gene_list))
length(gene_list)
ONT='BP'

gse <- clusterProfiler::gseGO(gene_list, 
                              ont=ONT, 
                              keyType = 'ENSEMBL', 
                              OrgDb = 'org.Hs.eg.db', 
                              pvalueCutoff  = 0.05)

#gse <- clusterProfiler::gseKEGG(gene_list, 
#                              ont=ONT, 
#                              keyType = 'ENSEMBL', 
#                              OrgDb = 'org.Hs.eg.db', 
#                              pvalueCutoff  = 0.05)
#
#
require(DOSE)
library('enrichplot')
#BiocManager::install('ggnewscale')
library('ggnewscale')


results_file<-paste0(outdir_s, '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T)
write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))

jpeg(paste0(results_file, '.jpeg'))
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()


x2 <- pairwise_termsim(gse)

p<-emapplot(x2,showCategory = 20,
            layout = "nicely")
p <- p + theme(text=element_text(size=12))
p

ggsave(paste0(results_file, '_emap',  '.jpeg'), width=8, height=8)





### 
library(ggplot2)
library(dplyr)
library(stringr)
gse@result
x=gse
## count the gene number
gene_count<- x@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)

## merge with the original dataframe
dot_df<- left_join(x@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

## plot
library(forcats) ## for reordering the factor
ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GO pathway enrichment")






gsea_run <- clusterProfiler::GSEA(gene_list, 
                              keyType = 'ALIAS', 
                              pAdjustMethod = 'BH')


