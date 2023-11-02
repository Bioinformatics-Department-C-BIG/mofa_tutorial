




process_mirnas=FALSE
source(paste0(script_dir, '/ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))
## 1. get se
## 2. deseq 
## 3. enrichment 
se_clusters<-se_filt_V08
deseq_all_groups <- vector("list", length = 3)


se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust )

se_clusters$kmeans_grouping
se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)

cd<-colData(se_clusters)
colData(se_clusters)[cd$INEXPAGE%in%'INEXHC','kmeans_grouping']<-'HC'
se_clusters$kmeans_grouping<-as.factor(se_clusters$kmeans_grouping)

se_filt

### TODO: add LEDD --> Turn it to zero for the ones that dont take it? 
se$LEDD


add_med='LEDD_scaled'
add_med='PDMEDYN'

if (add_med=='PDMEDYN'){
  formula_deseq = '~AGE_SCALED+SEX+PDMEDYN+kmeans_grouping'
  
}else if(add_med=='LEDD_scaled') {
  formula_deseq = '~AGE_SCALED+SEX+LEDD_scaled+kmeans_grouping'
  
}else{
  formula_deseq = '~AGE_SCALED+SEX+kmeans_grouping'
  
}
formula_deseq

deseq_by_group<-function(se_filt, formula_deseq){
  
        se_filt$kmeans_grouping
        se_filt<-preprocess_se_deseq2(se_filt)
        se_filt$kmeans_grouping
        se_filt$PDMEDYN = as.factor(se_filt$PDMEDYN)
        se_filt$PDMEDYN[is.na(se_filt$PDMEDYN)]=0
        
        
        se_filt$LEDD[is.na(se_filt$LEDD)]=0 # add zeros to na then scale!
        se_filt$LEDD_scaled<- scale(se_filt$LEDD)

        
    
        ddsSE <- DESeqDataSet(se_filt, 
                              design = as.formula(formula_deseq))
        ddsSE<-estimateSizeFactors(ddsSE)
        
        vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
        
        
        deseq2Data <- DESeq(ddsSE)
        deseq2Results<-results(deseq2Data)
        deseq2ResDF <- as.data.frame(deseq2Results)
        
        padj_T_hv<-0.05
        log2fol_T_hv<-0.1
        
        ### this is also done later on -- save from there? 
        deseq2ResDF$mofa_sign<- ifelse(deseq2ResDF$padj <padj_T_hv & abs(deseq2ResDF$log2FoldChange) >log2fol_T_hv , "Significant", NA)
        deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
        
        return(deseq2ResDF)
}

cluster_id=1
se_filt=se_clusters[,se_clusters$kmeans_grouping %in% c(cluster_id,'HC')]
deseq_all_groups[[cluster_id]]=deseq_by_group(se_filt, formula_deseq)

cluster_id=3
se_filt=se_clusters[,se_clusters$kmeans_grouping %in% c(cluster_id,'HC')]
deseq_all_groups[[cluster_id]]=deseq_by_group(se_filt, formula_deseq)





deseq2ResDF3<-deseq_all_groups[[3]]
deseq2ResDF1<-deseq_all_groups[[1]]
dim(deseq2ResDF3)
dim(deseq2ResDF1)

deseq_sig3<-deseq2ResDF3[deseq2ResDF3$mofa_sign %in% 'Significant',]
deseq_sig1<-deseq2ResDF1[deseq2ResDF1$mofa_sign %in% 'Significant',]


dim(deseq_sig1)
dim(deseq_sig3)


#union_de<-list(group1=de_group_vs_control1$symbol, group2=de_group_vs_control2$symbol, group3=de_group_vs_control3$symbol)
#fname_venn=paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/','venn_de_per_group_',keep_all_feats,'.png')







deseq2ResDF$log2FoldChange
order_by_metric='log2FoldChange'
order_by_metric='log2pval'

ONT='BP'
pvalueCutoff_sig=0.05



cluster_id=1

gene_list1<-get_ordered_gene_list(deseq2ResDF1,  order_by_metric, padj_T=1, log2fol_T=0 )
names(gene_list1)<-gsub('\\..*', '',names(gene_list1))
results_file_cluster=paste0(outdir, '/clustering/enrichment/gseGO',add_med, '_', ONT, '_', order_by_metric, 'clust', cluster_id)
gse1<-run_enrich_per_cluster(deseq2ResDF1, results_file_cluster,N_DOT=20, N_EMAP=30 )

cluster_id=3
results_file_cluster=paste0(outdir, '/clustering/enrichment/gseGO',add_med, '_', ONT, '_', order_by_metric, 'clust', cluster_id)
gse3<-run_enrich_per_cluster(deseq2ResDF3, results_file_cluster,N_DOT=20, N_EMAP=30)
gene_list3<-get_ordered_gene_list(deseq2ResDF3,  order_by_metric, padj_T=1, log2fol_T=0 )
names(gene_list3)<-gsub('\\..*', '',names(gene_list3))






gse_compare<-compareCluster(geneClusters = list(G1=gene_list1,G3=gene_list3 ), 
                            fun = "gseGO", 
                            OrgDb='org.Hs.eg.db', 
                            ont=ONT, 
                            keyType = 'ENSEMBL') 

### RUN SCRUPTI compare







