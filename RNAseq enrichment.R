
# BiocManager::install('org.Mm.eg.db')
#BiocManager::install('org.Hs.eg.db')

BiocManager::install('clusterProfiler', force=TRUE)
#BiocManager::install('apeglm')
# BiocManager::install('AnnotationDbi')

library(DESeq2)
#install.packages('contextual')
#BiocManager::install('DOSE', force=TRUE)


library(pheatmap)
#library(contextual)

library(org.Hs.eg.db)
library(DOSE)


library('apeglm')
 
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)

source('bladder_cancer/deseq2_vst_preprocessing.R')


isRStudio <- Sys.getenv("RSTUDIO") == "1"




#### Try using the deseq2 enrichment results 
### ADDD FILTERS HERE

res=deseq2ResDF

T_lfc=0.1
padj_T=0.05
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < padj_T , "Significant", NA)
deseq2ResDF$sign_lfc <- ifelse(deseq2ResDF$padj < padj_T & abs(deseq2ResDF$log2FoldChange) >T_lfc , "Significant", NA)


res=res[res$sign_lfc=='Significant'& !is.na(res$sign_lfc),]

dim(res)
res
# Order the DE gene list by the stat statistic 
#remove negatives thatw ere introduced with vst transofrmations
res<-res[res$baseMean>0,]
res<-res[res$baseMean>0,]

res <- res[order(-res$stat),]
gene_list<-res$stat
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

#ONT='MF'
#gse <- clusterProfiler::gseGO(gene_list, 
#                              ont=ONT, 
#                              keyType = 'ENSEMBL', 
#                              OrgDb = 'org.Hs.eg.db', 
#                              pvalueCutoff  = 0.05)
#
require(DOSE)

dev.off()

jpeg(paste0(outdir_s, '/gseGO', ONT, '_',padj_T, '_', T_lfc, '.jpeg'))
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
dev.off()


emapplot(gse, showCategory = 10)



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


