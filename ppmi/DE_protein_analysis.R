# Script to run limma for proteins 

source(paste0('ppmi/setup_os.R'))
library(limma)
library(pheatmap)
library(R.filesets)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(ggplot2)

process_mirnas=FALSE

VISIT='V06'
source(paste0(script_dir, 'ppmi/config.R' ))
source(paste0(script_dir,'ppmi/utils.R'))

source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))
source(paste0(script_dir,'/ppmi/plotting_utils.R'))


#BiocManager::install('rmdformats')
##### START HERE WITH PROTEOMICS 
## TODO: SAVE AND LOAD 
# se_filt and vsn mat 

### THIS IS ALREADY filtered by cohort and VISIT 
# 
dir.create(outdir_s_p)


datalist<-loadRDS(prot_vsn_se_filt_file)
vsn_mat<-datalist[[1]]
se_filt_prot<-datalist[[2]]

cohort_ids<-c(1,2)
names(cohort_ids)=c('INEXPD', 'INEXHC')

se_filt_prot$COHORT_orig<-se_filt_prot$COHORT


if (any(is.na(se_filt_prot$COHORT))){

  se_filt_prot$COHORT<- cohort_ids[se_filt_prot$INEXPAGE]

}
se_filt_prot$COHORT
tmp<- assays(se_filt_prot)[[1]]


if (run_vsn){
  protein_matrix_full<-vsn_mat
}else{
  protein_matrix_full<-tmp
  }
nfeats<-dim(protein_matrix_full)[1]

colnames(protein_matrix_full)

load_mofa_clusts=TRUE
MOFAobject_clusts
if (load_mofa_clusts){
    
        # MOFAobject_clusts<-MOFAobjectPD
        # Obtain clustering from mofa
        se_sel_proteins=se_filt_prot
        se_clusters_prot<-filter_se(se=se_sel_proteins, VISIT=VISIT, sel_coh = sel_coh, sel_sub_coh = sel_ps)
        se_clusters_prot$kmeans_grouping<- groups_from_mofa_factors(se_clusters_prot$PATNO, MOFAobject_clusts, y_clust );
        se_clusters_prot$kmeans_grouping=as.numeric(se_clusters_prot$kmeans_grouping)
        se_clusters_prot$kmeans_grouping[se_clusters_prot$COHORT==2]='HC'


}
se_clusters_prot$kmeans_grouping
clusters_indices
groups_from_mofa_factors(se_clusters_prot$PATNO, MOFAobject_clusts, y_clust )
  
  se_filt_all<-lapply(clusters_indices, function(cluster_id){ 
    se_return<-se_clusters_prot[,se_clusters_prot$kmeans_grouping %in% c(cluster_id)]
     if (length(unique(na.omit(se_return$COHORT)))==1){
         se_return$COHORT<-se_return$kmeans_grouping
     }
  
    return(se_return)

})
   


length(se_filt_all)

protein_matrices_all<-lapply(se_filt_all, function(se_cluster){
    protein_matrix_full  = protein_matrix_full[,colnames(protein_matrix_full) %in% se_cluster$PATNO_EVENT_ID ]
})


clust_id=3
results_de<-de_proteins_by_group(se_filt_all[[clust_id]], protein_matrices_all[[clust_id]])
protein_matrix = protein_matrices_all[[clust_id]]
se_filt<-se_filt_all[[clust_id]]; dim(se_filt); dim(protein_matrix)

results_de<-de_proteins_by_group(se_filt_all[[clust_id]], protein_matrices_all[[clust_id]])
 # write the de proteins
write.csv(results_de, paste0(outdir_s_p, 'results.csv'))







#### Enhanced volcano ####
library(EnhancedVolcano)
ylim=max(-log10(results_de$adj.P.Val))+0.5
x='logFC'
y= 'adj.P.Val'

prefix='prot_'
fname<-paste0(outdir_s_p,'/EnhancedVolcano_edited_', prefix,VISIT,'.jpeg')
pvol<-plotVolcano(results_de, se_filt, title='', xlim=NULL, lab=NULL,x=x, y= y,pCutoff = 10e-2,FCcutoff = 0.5)
pvol
ggsave(fname,pvol, width=6,height=8, dpi=300)




################### HEATMAPS  ############

#ARRANGE
order_by_metric<-'padj_reverse'
log2fol_T_overall=0;log2fol_T_hm=0.7;padj_T_overall=0.05;padj_T_hm=0.05 # cutoffs for the heatmaps 

# filter the gene list do add on the heatmap 
gene_list_limma_significant_hm=rownames(results_de)[results_de$adj.P.Val<padj_T_hm & 
                                                           abs(results_de$logFC)>log2fol_T_hm]


ids<-rownames(vsn_mat) %in% gene_list_limma_significant_hm

# heatmap input is the whole vsn_mat for ALL patients 
hm<-vsn_mat[ids,] # filter genes to load  




#### Create the heatmaps #### 

coldata_to_plot<-c('COHORT', 'SEX', 'AGE', 'PATNO_EVENT_ID', 'NP2PTOT' )
filter_highly_var=FALSE; most_var_t=FALSE
cluster_cols=TRUE
n_sig_f=50
fname<-paste0(outdir_s_p, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
              filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols, '.jpeg')

hm_plot<-create_heatmap_proteins(hm=hm,se_filt=se_filt, fname=fname, coldata_to_plot=coldata_to_plot)



# TODO: provide the cluster labels 



### Create a holistic map ####
# This map contains all top mofa genes, 
# all patients, and their cluster labels 
# top features to plot from mofa 

mofa_view=ifelse(TISSUE=='Plasma', 'proteomics_plasma', 'proteomics_csf')
# choose 11, 23 for proteins 
top_cluster_feats<-concatenate_top_features(MOFAobject_clusts, factors_all=c(11,23),view=mofa_view, top_fr=0.01 )
top_cluster_feats$feature <-gsub('\\_.*','',top_cluster_feats$feature);top_cluster_feats
top_cluster_feats
ids<-rownames(vsn_mat) %in% top_cluster_feats$feature

# heatmap input is the whole vsn_mat for ALL patients 
hm<-vsn_mat[ids,] # filter genes to load  

coldata_to_plot<-c('COHORT', 'SEX', 'AGE', 'PATNO_EVENT_ID', 'NP2PTOT', 'kmeans_grouping' )
filter_highly_var=FALSE; most_var_t=FALSE
cluster_cols=TRUE
n_sig_f=50
fname<-paste0(outdir_s_p, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
              filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,'all_clusters.jpeg')

#hm=hm;se_filt=se_clusters_prot
hm_plot<-create_heatmap_proteins(hm=hm,se_filt=se_clusters_prot, fname=fname, coldata_to_plot=coldata_to_plot)





################# ENRICHMENT - GSEA-GO #############
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'log2pval'
order_statistic<-'logFC'


log2fol_T_overall<-0.1
padj_T_overall<-.05



################### run gsea or ora with anova ######################
gene_list1<-results_de[,order_statistic]
names(gene_list1)<-rownames(results_de)
gene_list_ord=gene_list1[order(-gene_list1)]
gene_list_limma_significant=rownames(results_de)[results_de$adj.P.Val<padj_T_overall]
gene_list_limma_significant_pval=rownames(results_de)[results_de$P.Value<T]


run_ORA=TRUE; 
use_protein_pval=FALSE ## Proteins to use a sinput 
use_pval=TRUE  ### WHAT TO PLOT## these do not work well with false so keep it and mark the number of sig in the legend


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
length(gene_list_ora)
if (run_ORA){

    run_ora_gene_list(gene_list_ord, results_file_ora = outdir_s_p_enrich_file, keyType='SYMBOL', N_DOT=15, N_EMAP=25, N_NET=25,
    pvalueCutoff_sig = 0.05, top_p = 200)
}else{
    # TODO: run gsea 
    print('run gsea needs implementation')

}















