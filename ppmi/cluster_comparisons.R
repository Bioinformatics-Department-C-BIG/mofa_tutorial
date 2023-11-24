process_mirnas=FALSE
source(paste0(script_dir, '/ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))


## 1. get Summarized Experiment with metrics from all time points 
## 2. Run deseq 
## 3. enrichment 


MOFAobject_clusts=MOFAobjectPD


# TODO: make function to load for rnas and mirnas separately

se_clusters<-filter_se(se_filt_combat, VISIT='V08', sel_coh = sel_coh, sel_sub_coh = sel_ps) # se_filt_combat is missing one sample that was on one plate 
se_clusters<-filter_se(se_rnas, VISIT='V08', sel_coh = sel_coh, sel_sub_coh = sel_ps)
se_clusters<-filter_se(se_filt_V08, VISIT='V08', sel_coh = sel_coh, sel_sub_coh = sel_ps)


### Decide on the parameters settings 
# Set the outdirectory 

y_clust='NP2PTOT_LOG'
clust_name=paste0(y_clust, '_clust')

## Outputs 
# 1. DE files 
# 2. Venns 



formula_deseq = '~AGE_SCALED+SEX+kmeans_grouping'
#formula_deseq = '~AGE_SCALED+SEX+Plate+kmeans_grouping'

assay(se_clusters)


MOFAobject_clusts<-MOFAobjectPD
deseq_all_groups <- vector("list", length = 3);
se_filt_all<- vector("list", length = 3);
y_clust='NP2PTOT_LOG'
se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust )

se_clusters$kmeans_grouping
se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)


nclusts=length(table(se_clusters$kmeans_grouping));nclusts
cluster_params_dir<-paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/')

fname_venn=paste0(cluster_params_dir, '/','venn_de_per_group_deseq.png');fname_venn


cd<-colData(se_clusters)
colData(se_clusters)[cd$INEXPAGE%in%'INEXHC','kmeans_grouping']<-'HC'
se_clusters$kmeans_grouping<-as.factor(se_clusters$kmeans_grouping)

se_filt

### TODO: add LEDD --> Turn it to zero for the ones that dont take it? 
se$LEDD


add_med='LEDD_scaled'
add_med='PDMEDYN'
add_med=FALSE
if (add_med=='PDMEDYN'){
  formula_deseq = '~AGE_SCALED+SEX+PDMEDYN+kmeans_grouping'
  
}else if(add_med=='LEDD_scaled') {
  formula_deseq = '~AGE_SCALED+SEX+LEDD_scaled+kmeans_grouping'
  
}else{
  formula_deseq = '~AGE_SCALED+SEX+SITE+kmeans_grouping'
  formula_deseq = '~AGE_SCALED+SEX+Plate+kmeans_grouping'
  
  formula_deseq = '~AGE_SCALED+SEX+kmeans_grouping'
  
}

formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+Plate+kmeans_grouping'

deseq_by_group<-function(se_filt, formula_deseq){
  
  
  
        # TODO: add plate and/OR site 
        # se_filt1 neutrophil counts, and usable bases
        se_filt$SITE<-as.factor(se_filt$SITE)
        se_filt$Usable_Bases_SCALE<-scale(se_filt$`Usable.Bases....`)
        
        se_filt<-preprocess_se_deseq2(se_filt)
        se_filt$kmeans_grouping
        se_filt$PDMEDYN = as.factor(se_filt$PDMEDYN)
        se_filt$PDMEDYN[is.na(se_filt$PDMEDYN)]=0
        

        se_filt$LEDD[is.na(se_filt$LEDD)]=0 # add zeros to na then scale!
        se_filt$LEDD_scaled<- scale(se_filt$LEDD)

        ddsSE <- DESeqDataSet(se_filt, 
                              design = as.formula(formula_deseq))
        ddsSE<-estimateSizeFactors(ddsSE)
        
        #vsd <- varianceStabilizingTransformation(ddsSE, blind=FALSE)
        
        
        deseq2Data <- DESeq(ddsSE, parallel=TRUE)
        deseq2Results<-results(deseq2Data)
        deseq2ResDF <- as.data.frame(deseq2Results)
        
        padj_T_hv<-0.05
        log2fol_T_hv<-0.1
        
        ### this is also done later on -- save from there? 
        deseq2ResDF$mofa_sign<- ifelse(deseq2ResDF$padj <padj_T_hv & abs(deseq2ResDF$log2FoldChange) >log2fol_T_hv , "Significant", NA)
        deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
        

        return(deseq2ResDF)
}




deseq_all<- vector("list", length = 3)

### se_filt_all: list to hold the se
### deseq_all_groups: list to hold the deseq results 
### deseq_significant_all_groups: list to hold significant 


for (cluster_id in 1:3){
  

  ### 1. for each cluster, create se filt with controls, 
  ### 2. run deseq 
  ### 3. get significant per cluster 

  de_file<-paste0(cluster_params_dir, '/de_cluster_', cluster_id , '.csv')

  se_filt_all[[cluster_id]]<-se_clusters[,se_clusters$kmeans_grouping %in% c(cluster_id,'HC')]

  # if deseq exists load:
  if (file.exists(de_file)){
    # if de file exists load it - unfiltered de results file
    deseq2ResDF<-read.csv(paste0(de_files,  cluster_id , '.csv'), row.names=1 )

  }else{
    # else run the deseq with the design formula specified 
        deseq2ResDF=deseq_by_group(se_filt_all[[cluster_id]], formula_deseq)
        deseq_all_groups[[cluster_id]]<-deseq2ResDF
        deseq2ResDF$GENE_SYMBOL<-get_symbols_vector(gsub('\\..*', '',rownames(deseq2ResDF)))
        write.csv(deseq2ResDF, de_files[[cluster_id]], row.names=TRUE)
      
  }
  deseq_all_groups[[cluster_id]]<-deseq2ResDF
  deseq_all[[cluster_id]]<-deseq2ResDF[deseq2ResDF$mofa_sign %in% 'Significant',] # holds the significant only

} 

# Save and load 

deseq_all_names<-lapply(deseq_all, function(x){return(  gsub('\\..*', '',rownames(x))   )  })
names(deseq_all_names)<-paste0('SG', 1:length(deseq_all_names))




#### 1. Venn from significant 
fname_venn
create_venn(venn_list = deseq_all_names, fname_venn =fname_venn,main =paste0( ' DE molecules for each molecular cluster' ))

graphics.off()
# TODO: 




########### 



# TODO: venn before and after correction 
plate

cluster_id = 1

se_filt=se_filt_all[[cluster_id]]

deseq2ResDF=deseq_all_groups[[cluster_id]]

#any(rownames(deseq2ResDF) %in% remove_genes)
pvol<-plotVolcano(deseq2ResDF, se_filt, title=paste0('Cluster ', cluster_id), xlim=c(-1.1,1.1))

outdir_s
pvol
fname
fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT,'.jpeg')
fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT_S, '_cluster_',cluster_id, '.jpeg')

ggsave(fname,pvol, width=9,height=12, dpi=300)








#### 2. Venn from significant in top of factor

## 1. get top of factor / or all higly variable genes input into MOFA
## 2. intersect with DE 
sel_factor=2
top10<- gsub('\\..*', '',select_top_bottom_perc(MOFAobject=MOFAobject, view=2, factors = sel_factor, top_fr = 0.1))
top10<- gsub('\\..*', '',select_top_bottom_perc(MOFAobject=MOFAobject, view=2, factors = sel_factor, top_fr = 1))
length(top10)
# intersect with the top factors 
deseq_all_top<-lapply(deseq_all_names, function(x) intersect(x,top10 ) )


deseq_all_top


fname_venn=paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/','venn_de_per_group_deseq', 'top_f', sel_factor ,  '.png')
create_venn(venn_list = deseq_all_top, fname_venn =fname_venn,main =paste0( ' DE molecules for each molecular cluster AND highly variable' ))




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







