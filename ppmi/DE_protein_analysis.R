# Script to run limma for proteins 

source(paste0('ppmi/setup_os.R'))
library(limma)
library(pheatmap)
library(R.filesets)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)

process_mirnas=FALSE

VISIT='V06'
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
prot_vsn_se_filt_file
datalist<-loadRDS(prot_vsn_se_filt_file)
vsn_mat<-datalist[[1]]
se_filt<-datalist[[2]]

cohort_ids<-c(1,2)
names(cohort_ids)=c('INEXPD', 'INEXHC')

se_filt$COHORT_orig<-se_filt$COHORT


if (any(is.na(se_filt$COHORT))){

  se_filt$COHORT<- cohort_ids[se_filt$INEXPAGE]

}
se_filt$COHORT
tmp<- assays(se_filt)[[1]]


if (run_vsn){
  protein_matrix<-vsn_mat
}else{
  protein_matrix<-tmp
  }
nfeats<-dim(protein_matrix)[1]


protein_matrix




results_de<-de_proteins_by_group(se_filt, protein_matrix)
#results_de<-topTable(fit.cont, coef='COHORT' )

#FC= mean(condition_A_replicates)n/ mean(control_replicates)n   




dir.create(outdir_s_p)

ns_full<-table(se_filt$COHORT_DEFINITION)
ns<-paste0(rownames(ns_full)[1],' ', ns_full[1], '\n' ,names(ns_full)[2], ' ', ns_full[2])

# TODO: enhanced volcano
library(EnhancedVolcano)
ylim=max(-log10(results_de$adj.P.Val))+0.5
ylim
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
fname<-paste0(outdir_s_p,'/EnhancedVolcano_edited_', prefix,VISIT,'.jpeg')
ggsave(fname,pvol, width=6,height=8, dpi=300)

#ggsave(fname,pvol, width=6,height=8, dpi=300)


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
if (TISSUE=='CSF' & VISIT=='V08'){
  log2fol_T_overall=1
  
}
log2fol_T_overall=0
log2fol_T_hm=0.7
padj_T_overall=0.05
padj_T_hm=0.05

results_de
results_de$logFC
gene_list_limma_significant_heatmap=rownames(results_de)[results_de$adj.P.Val<padj_T_hm & 
                                                           abs(results_de$logFC)>log2fol_T_hm]
gene_list_limma_significant_heatmap

length(gene_list_limma_significant_heatmap)
ids<-rownames(vsn_mat) %in% gene_list_limma_significant_heatmap


hm<-vsn_mat[ids,]
results_de$padj_reverse<--results_de$adj.P.Val




df<-as.data.frame(colData(se_filt)[c('COHORT', 'SEX', 'AGE', 'PATNO_EVENT_ID' )]); rownames(df)<-df$PATNO_EVENT_ID
df$PATNO_EVENT_ID<-NULL

#hm_ord<-hm[,order(df$COHORT)]

#<-paste0(outdir_s_p, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, '_', n_sig_f,'.jpeg')
filter_highly_var=FALSE; most_var_t=FALSE
cluster_cols=TRUE
n_sig_f=50
fname<-paste0(outdir_s_p, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
              filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols, '.jpeg')
fname
graphics.off()
library(ggplot2)
jpeg(fname, res = 200, width=3+log2(dim(hm)[2]), height=4+log2(dim(hm)[1]), units='in')

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
fname






################# ENRICHMENT - GSEA-GO #############
order_statistic<-'log2pval'
order_statistic<-'logFC'
order_statistic<-'log2pval'
order_statistic<-'logFC'



#order_statistic<-'log2pval_not_adj' - NO RESULTS 
results_de$pval_reverse<- -results_de$P.Value ## this prob does not work with gsea because there are no negative values 

results_de$log2pval<- -log10(results_de$adj.P.Val) * results_de$logFC
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


highly_variable_sign_proteins_mofa<-highly_variable_proteins_mofa[rownames(highly_variable_proteins_mofa) %in%  signif_proteins_mofa_ids,]


results_de_signif$abslog2pval

################### run gsea or ora with anova ######################
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
}












