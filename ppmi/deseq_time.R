

###
#BiocManager::install('DEGreport')
install.packages('lasso2')
#library('lasso2')
library('DEGreport')

# https://github.com/hbctraining/DGE_workshop_salmon_online/blob/master/lessons/08a_DGE_LRT_results.md

### aply deseq time PER GROUP
resultsNames(deseq2Data)


res_LRT <- results(deseq2Data)
res_LRT$padj
deseq2Results_time <- results(deseq2Data,  contrast=c('EVENT_ID', 'V08','BL'))
deseq2Results_time <- results(deseq2Data,  contrast=c('EVENT_ID', 'V08','BL'))


deseq2Data$COHORT <- relevel(deseq2Data$COHORT, "1")
resultsNames(deseq2Data)

deseq2Results_time <- results(deseq2Data)

deseq2Results_time
#deseq2Results_time <- results(deseq2Data,  contrast=c('COHORT', '1','2'))

deseq2Results_time$padj



# Create a tibble for LRT results
res_LRT_tb <- deseq2Results_time %>%
  data.frame() %>%
  tibble::rownames_to_column(var='gene') %>%
  as_tibble()


res_LRT_tb$padj
# Subset to return genes with padj < 0.05
padj.cutoff=0.00001
padj.cutoff=0.05

sigLRT_genes <- res_LRT_tb %>% 
  dplyr::filter(padj < padj.cutoff)

# Get number of significant genes
nrow(sigLRT_genes)
sigLRT_genes

# Compare to numbers we had from Wald test


clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)

cluster_rlog<-vsd[sigLRT_genes$gene, ]

clusters <- degPatterns(cluster_rlog, metadata = meta, time = "EVENT_ID", col=NULL)
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "COHORT", col=NULL)



#### TODO: REMOVE OTHER time points and also try deseq with patients ids
