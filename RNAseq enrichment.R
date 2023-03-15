
# BiocManager::install('org.Mm.eg.db')
 BiocManager::install('org.Hs.eg.db')

# BiocManager::install('clusterProfiler')
# BiocManager::install('apeglm')
 BiocManager::install('AnnotationDbi')

library(DESeq2)
library(pheatmap)

library(org.Hs.eg.db)
library(DOSE)
library(pathview)


library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)
library(tidyverse)

source('bladder_cancer/deseq2_vst_preprocessing.R')
condition=sample_info$Subtype
dds_enrich = DESeqDataSet(dds, design=~Subtype)

# This function runs differential exppression with subtype as the two groups 

dds_enrich = DESeq(dds_enrich)

resultsNames(dds_enrich)


res <- results(dds_enrich, name = 'Subtype_NPS3_vs_NPS1')
resLFC <- lfcShrink(dds_enrich, coef="Subtype_NPS3_vs_NPS1", type="apeglm")
res

res = results(dds_enrich, contrast = c('Subtype', 'NPS1', 'NPS3' ))

# Order the DE gene list by the stat statistic 
#remove negatives thatw ere introduced with vst transofrmations
res<-res[res$baseMean>0,]

res <- res[order(-res$stat),]

gene_list<-res$stat
names(gene_list)<-rownames(res)
gene_list


top20<-res[1:200,]
write.csv(top20, 'bladder_cancer/Enrichment/top20.txt')
#  Takes input the DE genes from DESeq2 
gse <- clusterProfiler::gseGO(gene_list, 
             ont='BP', 
             keyType = 'ALIAS', 
             OrgDb = 'org.Hs.eg.db')


require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
emapplot(gse, showCategory = 10)



gsea_run <- clusterProfiler::GSEA(gene_list, 
                              keyType = 'ALIAS', 
                              pAdjustMethod = 'BH')


