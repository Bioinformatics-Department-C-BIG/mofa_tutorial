
library(DESeq2)
library(pheatmap)

library(org.Hs.eg.db)
#BiocManager::install('clusterProfiler')
library(DOSE)
require(DOSE)


library('apeglm')
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)
library(tidyverse)

#source('bladder_cancer/deseq2_vst_preprocessing.R')
#condition=sample_info$Subtype
#dds_enrich = DESeqDataSet(dds, design=~Subtype)

# This function runs differential exppression with subtype as the two groups 

#dds_enrich = DESeq(dds_enrich)

#resultsNames(dds_enrich)


#res <- results(dds_enrich, name = 'Subtype_NPS3_vs_NPS1')
#resLFC <- lfcShrink(dds_enrich, coef="Subtype_NPS3_vs_NPS1", type="apeglm")

#res = results(dds_enrich, contrast = c('Subtype', 'NPS1', 'NPS3' ))

# Order the DE gene list by the stat statistic 
#remove negatives that were introduced with vst transofrmations




#### Try using the deseq2 enrichment results 

padj_T=0.05;log2fol_T=0.10
res=deseq2ResDF
res$sign_lfc <- ifelse(res$padj <padj_T & abs(res$log2FoldChange) >log2fol_T , "Significant", NA)

length(which(!is.na(res$sign_lfc )))
res=res[res$sign_lfc=='Significant'& !is.na(res$sign_lfc),]
dim(res)
res
# Order the DE gene list by the stat statistic 
#remove negatives thatw ere introduced with vst transofrmations


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



results_file<-paste0(outdir_s, '/gseGO', '_', ONT, '_', padj_T, '_',  log2fol_T)
write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))
jpeg(paste0(results_file, '.jpeg'))
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
dev.off()


emapplot(gse, showCategory = 10)



gsea_run <- clusterProfiler::GSEA(gene_list, 
                              keyType = 'ALIAS', 
                              pAdjustMethod = 'BH')


