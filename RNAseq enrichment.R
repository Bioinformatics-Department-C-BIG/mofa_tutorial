
# BiocManager::install('org.Mm.eg.db')
#BiocManager::install('org.Hs.eg.db')

#BiocManager::install('clusterProfiler')
#BiocManager::install('apeglm')
# BiocManager::install('AnnotationDbi')

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



#### Try using the deseq2 enrichment results 


res=deseq2ResDF
res=res[res$sign_lfc=='Significant'& !is.na(res$sign_lfc),]
dim(res)
res
resLFC# Order the DE gene list by the stat statistic 
#remove negatives thatw ere introduced with vst transofrmations
res<-res[res$baseMean>0,]
res<-res[res$baseMean>0,]

res <- res[order(-res$stat),]
gene_list<-res$stat
names(gene_list)<-rownames(res)


#  Takes input the DE genes from DESeq2



names(gene_list)<-gsub('\\..*', '',names(gene_list))
length(gene_list)

gse <- clusterProfiler::gseGO(gene_list, 
                              ont='BP', 
                              keyType = 'ENSEMBL', 
                              OrgDb = 'org.Hs.eg.db', 
                              pvalueCutoff  = 0.05)

ONT='MF'
gse <- clusterProfiler::gseGO(gene_list, 
                              ont=ONT, 
                              keyType = 'ENSEMBL', 
                              OrgDb = 'org.Hs.eg.db', 
                              pvalueCutoff  = 0.05)

require(DOSE)
gse
dev.off()

jpeg(paste0(outdir_s, '/gseGO', ONT, '.jpeg'))
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()


emapplot(gse, showCategory = 10)



gsea_run <- clusterProfiler::GSEA(gene_list, 
                              keyType = 'ALIAS', 
                              pAdjustMethod = 'BH')


