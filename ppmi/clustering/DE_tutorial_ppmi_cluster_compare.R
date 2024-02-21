# Script to run limma for proteins 
# TODO: add venn per cluster 
# TODO: run pathways per cluster 

source(paste0('ppmi/setup_os.R'))
library(limma)
library(pheatmap)
library(R.filesets)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)

process_mirnas=FALSE

#VISIT_COMP='V06' # set elsewhere
VISIT = VISIT_COMP
source(paste0(script_dir, 'ppmi/config.R' ))
source(paste0(script_dir,'ppmi/utils.R'))

source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))


#BiocManager::install('rmdformats')
##### START HERE WITH PROTEOMICS 
## TODO: SAVE AND LOAD 
# se_filt and vsn mat 

### THIS IS ALREADY filtered by cohort and VISIT 
# 

prot_vsn_se_filt_file
datalist<-loadRDS(prot_vsn_se_filt_file)
vsn_mat<-datalist[[1]]
se_filt_proteins<-datalist[[2]]

cohort_ids<-c(1,2)
names(cohort_ids)=c('INEXPD', 'INEXHC')

se_filt_proteins$COHORT_orig<-se_filt_proteins$COHORT

# fix the cohort
if (any(is.na(se_filt_proteins$COHORT))){

  se_filt_proteins$COHORT<- cohort_ids[se_filt_proteins$INEXPAGE]

}
tmp<- assays(se_filt_proteins)[[1]]


if (run_vsn){
  protein_matrix_full<-vsn_mat
}else{
  protein_matrix_full<-tmp
  }
nfeats<-dim(protein_matrix_full)[1]

se_clusters<-se_filt_proteins


### attach mofa clusts to se
# MOFAobject_clusts<-MOFAobjectPD
# Obtain clustering from mofa


sm<-samples_metadata(MOFAobject_clusts)

sm$NP2PTOT_LOG_clust
patnos=se_clusters$PATNO
se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust );
#se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)
se_clusters$kmeans_grouping

protein_matrices<-list()

de_all_groups_proteins<-list()


# TODO: ensure that this is also run for all subjects together

length(se_filt_all)

for (cluster_id in 1:3){

  se_clusters
  se_cluster_ind<-se_clusters$kmeans_grouping %in% c(cluster_id,'HC')
  se_filt_all[[cluster_id]]<-se_clusters[,se_cluster_ind]
  protein_matrices[[cluster_id]]<-protein_matrix_full[,se_cluster_ind ]

}

protein_matrix<-protein_matrices[[1]]
se_filt_all[[1]]$COHORT



for (cluster_id in 1:3){

        outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )
        outdir_s_p
        se_filt_clust<-se_filt_all[[cluster_id]]
        protein_matrix<-protein_matrices[[cluster_id]]
        results_de<-de_proteins_by_group(se_filt=se_filt_clust, protein_matrix=protein_matrix)
        #results_de<-topTable(fit.cont, coef='COHORT' ) 
        de_all_groups_proteins[[cluster_id]]<-results_de
        #FC= mean(condition_A_replicates)n/ mean(control_replicates)n   


    #print(results_de['IL5',])
    #print(results_de['MMP9',])


    dir.create(outdir_s_p)

    ns_full<-table(se_filt_clust$COHORT_DEFINITION)
    ns<-paste0(rownames(ns_full)[1],' ', ns_full[1], '\n' ,names(ns_full)[2], ' ', ns_full[2])

    # TODO: enhanced volcano
    library(EnhancedVolcano)
    ylim=max(-log10(results_de$adj.P.Val))+0.5
    ylim
    #results_de<-de_all_groups_proteins[[3]]
    pvol<-EnhancedVolcano(results_de, 
                    lab = rownames(results_de),
                    x = 'logFC',
                    y = 'adj.P.Val', 
                    pCutoff = 10e-2,
                    FCcutoff = 0.5, 
                    ylim=c(0,ylim), 
                    title='', 
                    subtitle=ns
    )
    pvol
    prefix='prot_'
    fname_vol<-paste0(outdir_s_p,'/Volcano_', prefix, TISSUE,'_',VISIT_COMP,'_cluster_',cluster_id, '.jpeg')
    ggsave(fname_vol,pvol, width=6,height=8, dpi=300)
}

results_de
#ggsave(fname,pvol, width=6,height=8, dpi=300)

print(fname_vol)
## Create a p-adjusted
## how many total proteins?
dim(vsn_mat)


dim(vsn_mat)[1]

#common_de<-intersect(all_sig_proteins,anova_results_oneway_significant)
#dir.create(outdir_s_p)
outdir_s_p_enrich<-paste0(outdir_s_p, '/enrichment/'); dir.create(outdir_s_p_enrich)
#write.csv(common_de, paste0(outdir_s_p, 'common_de.csv'))



#fit.cont_sig[common_de]
#



################### HEATMAPS  ############

#ARRANGE
#df_ord<-df[order(df$COHORT),]
pvol
order_by_metric<-'padj_reverse'
if (TISSUE=='CSF' & VISIT_COMP=='V08'){
  log2fol_T_overall=1
  
}
log2fol_T_overall=0
log2fol_T_hm=0
padj_T_overall=0.05
padj_T_hm=0.05

results_de
gene_list_limma_significant_heatmap=rownames(results_de)[results_de$adj.P.Val<padj_T_hm & 
                                                           results_de$logFC>log2fol_T_hm]


length(gene_list_limma_significant_heatmap)
ids<-rownames(vsn_mat) %in% gene_list_limma_significant_heatmap
## if zero significant plot the highly variable proteins 
highly_variable_proteins_mofa<-  selectMostVariable(vsn_mat, TOP_PN)
ids_highly_var<-ids<-rownames(vsn_mat) %in% rownames(highly_variable_proteins_mofa)
ids_highly_var

hm<-vsn_mat[ids,]
if (sum(ids)<1){
  hm<-vsn_mat[ids_highly_var,]
}

results_de$padj_reverse<--results_de$adj.P.Val



# se_filt proteins contains them all?
se_filt_proteins$PATNO_EVENT_ID
df<-as.data.frame(colData(se_filt_proteins)[c('COHORT', 'SEX', 'AGE', 'PATNO_EVENT_ID' )]); rownames(df)<-df$PATNO_EVENT_ID
df$PATNO_EVENT_ID<-NULL
colnames(vsn_mat)
rownames(df)
#rownames(df)<-se_filt$PATNO_EVENT_ID
se_filt_proteins$COHORT
dim(df)
dim(hm)


#<-paste0(outdir_s_p, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, '_', n_sig_f,'.jpeg')
filter_highly_var=FALSE; most_var_t=FALSE


cluster_cols=TRUE
n_sig_f=30
fname<-paste0(outdir_s_p, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
              filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols, '.jpeg')

graphics.off()
library(ggplot2)

#if(process_mirnas){
  #lab=rownames(rowData(vsd_filt_genes)) }else{
   # lab=as.character(rowData(vsd_filt_genes)$SYMBOL)}

        #jpeg(fname, width=10*100, height=10*100, res=150)


        my_pheatmap<-pheatmap(hm, 
                              #labels_row=lab,
                              cluster_rows=TRUE, 
                              show_rownames=TRUE,
                              scale='row', 
                              cluster_cols=cluster_cols,
                              annotation_col=df, 
                              clustering_method = 'complete', 
                              clustering_distance_cols = 'euclidean'
        )
        dim(hm)

show(my_pheatmap)
dev.off()
my_pheatmap
fname
#ggsave(fname, my_pheatmap, dpi = 200, width=dim(hm)[2]/5+2, height=dim(hm)[1]/10+4)






#############################
#install.packages("OlinkAnalyze")
run_olink=FALSE
if (run_olink){
  
    
    library(OlinkAnalyze)
    #### ALSO TRY T-TEST
    data<-read_NPX(prot_files[1])
    data=prot_bl
    VISIT_COMP='BL'
    colnames(data)<- c( "PATNO" ,       "EVENT_ID"  ,   "Index"   ,     "OlinkID"    ,
                        "UniProt"     , "Assay"     ,   "MISSINGFREQ", "Panel"    , 
                        "PANEL_LOT_NR" ,"PLATEID"   ,   "QC_WARNING" ,  "LOD"    ,    
                        "NPX"   ,       "update_stamp")
    
    table(data$EVENT_ID)
    data_filt<-data[data$EVENT_ID %in% VISIT_COMP,]
    data_filt$EVENT_ID
    
    unique(data_filt$PATNO)
    outcome<-combined[,c('PATNO','COHORT', 'AGE_SCALED', 'SEX' )]
    outcome<-outcome[!duplicated(outcome$PATNO),]
    outcome<-outcome[outcome$COHORT %in% sel_coh,]
    
    ### now merge the specific data visit and metadata
    data_merged<-merge(data_filt, outcome, by='PATNO')
    data_merged$COHORT<-as.factor(data_merged$COHORT)
    
    data_merged$SampleID<-data_merged$PATNO
    ttest_results<-olink_ttest(df = data_merged,
                variable = 'COHORT')
    
    olink_wilcox_results<-olink_wilcox(df = data_merged,
                 variable = 'COHORT')
    
    olink_wilcox_results%>%
      dplyr::filter(Threshold == 'Significant')
    
    
    library(dplyr)
    data_merged$SEX=as.factor(data_merged$SEX)
    data_merged_rescaled<-data_merged
    data_merged_rescaled$NPX=data_merged_rescaled$NPX*1e6
    anova_results_oneway <- olink_anova(df = data_merged_rescaled, 
                                        variable = c('COHORT' ,'AGE_SCALED', 'SEX'))
    
    data_merged$COHORT=as.factor(data_merged$COHORT)
    anova_results_oneway <- olink_anova(df = data_merged, 
                                        variable = 'COHORT',
                                        covariates = c('AGE_SCALED', 'SEX'))
    
    anova_results_oneway_significant <- anova_results_oneway %>%
      dplyr::filter(Threshold == 'Significant')
    
    
    anova_results_oneway_significant_genes <- anova_results_oneway %>%
      dplyr::filter(Threshold == 'Significant') %>%
      dplyr::pull(Assay, Adjusted_pval)
    
    length(anova_results_oneway_significant)
    
    anova_results_oneway_significant[1:20]
    anova_results_oneway$Assay
    write.csv(anova_results_oneway,paste0(output_files, 'olink_de.csv' ))
    
    
    anova_posthoc_oneway_results <- olink_anova_posthoc(df = data_merged, 
                                        variable = 'COHORT',
                                        covariates = c('AGE_SCALED', 'SEX'))
    
    


}

T=0.05



################# ENRICHMENT - GSEA-GO #############
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'log2pval'
order_statistic<-'adj.P.Val'
order_statistic<-'signlog2pval'
order_statistic<-'P.Value'
order_statistic<-'log2pval'
#order_statistic<-'pval_reverse'
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'logFC'
order_statistic<-'signlog2pval'
order_statistic<-'logFC'



#order_statistic<-'log2pval_not_adj' - NO RESULTS 
results_de$pval_reverse<- -results_de$P.Value ## this prob does not work with gsea because there are no negative values 

results_de$log2pval<- -log10(results_de$adj.P.Val) * results_de$logFC
results_de$abslog2pval<- abs(results_de$log2pval)
results_de$abslog2pval<- abs(results_de$log2pval)
results_de$signlogFCpval<- abs(results_de$log2pval)

results_de$log2pval_not_adj<- -log10(results_de$P.Value) * results_de$logFC
results_de$signlog2pval<- -log10(results_de$adj.P.Val) * sign(results_de$logFC)
results_de[results_de$adj.P.Val<0.05,'signlog2pval']

write.csv(results_de, paste0(outdir_s_p, 'results.csv'))


log2fol_T_overall<-0.1
padj_T_overall<-.05
results_de_signif<-mark_significant(results_de,padj_T = padj_T_overall, log2fol_T = log2fol_T_overall, 
                            padj_name ='adj.P.Val',log2fc_name = 'logFC' , outdir_single = outdir_s_p  )


log2fol_T_mofa<-0
padj_T_mofa<-.05
signif_proteins_mofa<-mark_significant(results_de,padj_T = padj_T_mofa, log2fol_T = log2fol_T_mofa, 
                                    padj_name ='adj.P.Val',log2fc_name = 'logFC' , outdir_single = outdir_s_p  )

signif_proteins_mofa_ids<-rownames(signif_proteins_mofa %>%
                            dplyr::filter(significant =='Significant'))




results_de_signif$abslog2pval

################### run gsea with anova ######################
gene_list1<-results_de[,order_statistic]
names(gene_list1)<-rownames(results_de)
gene_list_ord=gene_list1[order(-gene_list1)]
gene_list_ord
gene_list_limma_significant=rownames(results_de)[results_de$adj.P.Val<padj_T_overall]
gene_list_limma_significant_pval=rownames(results_de)[results_de$P.Value<T]


run_anova=FALSE
run_ORA=TRUE; 
use_protein_pval=FALSE ## Proteins to use a sinput 
use_pval=TRUE  ### WHAT TO PLOT## these do not work well with false so keep it and mark the number of sig in the legend


pvalueCutoff_sig=0.1
if (run_anova){
  order_statistic<-'statistic'
  gene_list2<-pull(anova_results_oneway, statistic)
  names(gene_list2)<-anova_results_oneway$Assay
  gene_list_ord=gene_list2[order(-gene_list2)]
  
  
}
length(gene_list_ord)
gene_list_ora=gene_list_limma_significant
if (use_protein_pval){
  gene_list_ora=gene_list_limma_significant_pval
  
}
gene_list_ora
ONT='BP'

gene_list_ord

pvalueCutoff=1
#outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT, '_', order_statistic)



outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT,  '_', order_statistic, '_ora_', run_ORA, 'ppval_', use_protein_pval, 
                               '_anova_', run_anova, 'pval_', use_pval )
res_path<-paste0(outdir_s_p_enrich_file, 'gse.RDS')


#if (file.exists(res_path)){
#  gse_protein_full=loadRDS(res_path)
#}else{

# TODO: run per cluster 
# Input:  gene_list_ord per cluster (coming from DE limma results )
# 
        if (run_ORA){
          # writes 
           gse_protein_full = run_ora_gene_list(  gene_list_ord = gene_list_ord, results_file_ora =outdir_s_p_enrich_file  )
          
        }else{
          
        
          
          #### TODO: ALSO RUN ENRICHMENT using ANOVA from 
          gse_protein_full <- clusterProfiler::gseGO(gene_list_ord, 
                                        ont=ONT, 
                                        keyType = 'SYMBOL', 
                                        OrgDb = 'org.Hs.eg.db', 
                                        pvalueCutoff  = pvalueCutoff 
                                          )
      
          
          
          use_protein_pval=FALSE # if we are rrunning gsea this is not used actually so set to false everytime
          
        }
  saveRDS(gse_protein_full, res_path)
  
  
#}


