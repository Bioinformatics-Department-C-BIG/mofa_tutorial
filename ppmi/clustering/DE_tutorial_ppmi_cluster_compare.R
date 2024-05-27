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
#VISIT_COMP = 'V08'

source(paste0(script_dir, 'ppmi/config.R' ))
source(paste0(script_dir,'ppmi/utils.R'))

source(paste0(script_dir,'/bladder_cancer/preprocessing.R'))


#BiocManager::install('rmdformats')
##### START HERE WITH PROTEOMICS 
## TODO: SAVE AND LOAD 
# se_filt and vsn mat 

### THIS IS ALREADY filtered by cohort and VISIT 
# 
# todo: CREATE SE FILT FILES FOR PROTEINS 
 # olink or untargeted 
 # TODO: fix some CLUSTERS have zero patients 

if (prot_de_mode=='t'){
    datalist<-loadRDS(prot_vsn_se_filt_file)
    vsn_mat<-datalist[[1]]
    se_filt_proteins<-datalist[[2]]
    tissue = TISSUE
}else{
   # proteins_un_plasma<-as.matrix(read.csv2(prot_untargeted_plasma_vsn_f,row.names=1, header=TRUE, check.names = FALSE))
    prot_un_vsn<-as.matrix(read.csv2(prot_untargeted_un_vsn_f, row.names=1, header=TRUE, check.names = FALSE))
    tissue = tissue_un
    vsn_mat = prot_un_vsn
    prot_vsn_se<-getSummarizedExperimentFromAllVisits(vsn_mat, combined_bl_log) # create sumarized experiment 
    se_filt_prot_vsn<- filter_se(prot_vsn_se, VISIT, sel_coh)
    vsn_mat<-assay(se_filt_prot_vsn) # reproduce after filtering by visit
    se_filt_proteins <- se_filt_prot_vsn

}

cohort_ids<-c(1,2)
names(cohort_ids)=c('INEXPD', 'INEXHC')

se_filt_proteins$COHORT_orig<-se_filt_proteins$COHORT
se_filt_proteins$INEXPAGE
# fix the cohort
if (any(is.na(se_filt_proteins$COHORT))){
  se_filt_proteins$COHORT<- cohort_ids[se_filt_proteins$INEXPAGE]
}
tmp<- assays(se_filt_proteins)[[1]]

head(as.matrix(vsn_mat));
head(as.matrix(tmp))
if (run_vsn){
  protein_matrix_full<-vsn_mat
}else{
  protein_matrix_full<-tmp
  }
nfeats<-dim(protein_matrix_full)[1]

se_clusters<-se_filt_proteins
se_filt_proteins$PATNO

### attach mofa clusts to se
# MOFAobject_clusts<-MOFAobjectPD
# Obtain clustering from mofa


sm<-samples_metadata(MOFAobject_clusts)

sm$NP2PTOT_LOG_clust
patnos=se_clusters$PATNO
y_clust<-DIFF_VAR

fact <- get_factors_for_metric(DIFF_VAR)
cluster_params_dir<-get_cluster_params_dir(DIFF_VAR)

se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust );
#se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)
se_clusters$kmeans_grouping

protein_matrices<-list()

de_all_groups_proteins<-list()


# TODO: ensure that this is also run for all subjects together

length(se_filt_all)
clusters_indices= list( c('1', 'HC'),c('2', 'HC'),c('3', 'HC'), c('1','2','3', 'HC')) 




for (cluster_id_num in 1:length(clusters_indices)){



  print(paste('cluster:',cluster_id_num))
  cluster_id_index=clusters_indices[[cluster_id_num]]
  #cluster_id_name=clusters_names[[cluster_id_num]]


  se_cluster_ind<-se_clusters$kmeans_grouping %in% c(cluster_id_index,'HC')

  
  se_filt_all[[cluster_id_num]]<-se_clusters[,se_cluster_ind]
  se_filt_all[[1]]$COHORT
  protein_matrices[[cluster_id_num]]<-protein_matrix_full[,se_cluster_ind ]

}

protein_matrix<-protein_matrices[[1]]
se_filt_all[[1]]$COHORT

jpeg(paste0(outdir, '/boxplots.jpeg'))
boxplot(protein_matrix)
dev.off()
prefix='prot_'
prefix_full =  paste0(prefix, tissue,'_', prot_de_mode)
print(cluster_params_dir)


for (cluster_id_num in  c(1:4)){
      print(cluster_id_num)
        cluster_id_index=clusters_indices[[cluster_id_num]]

        outdir_s_p <- paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' )

        se_filt_clust<-se_filt_all[[cluster_id_num]]
        protein_matrix<-protein_matrices[[cluster_id_num]]
        results_de<-de_proteins_by_group(se_filt=se_filt_clust, protein_matrix=protein_matrix)
        #results_de<-topTable(fit.cont, coef='COHORT' ) 
        de_all_groups_proteins[[cluster_id_num]]<-results_de

        #FC= mean(condition_A_replicates)n/ mean(control_replicates)n   




    suppressWarnings(dir.create(outdir_s_p, recursive=TRUE))
   rfile<-paste0(outdir_s_p, prefix_full, '_de_cl',cluster_id_num,  '_results.csv')
   print(rfile)
    write.csv(results_de, rfile)

    ns_full<-table(se_filt_clust$COHORT_DEFINITION)
    ns<-paste0(rownames(ns_full)[1],' ', ns_full[1], '\n' ,names(ns_full)[2], ' ', ns_full[2])

    # TODO: enhanced volcano
    library(EnhancedVolcano)
    ylim=max(-log10(results_de$adj.P.Val))+0.5
    xlim=max(abs(results_de$logFC))+0.1
    colnames(results_de)



   if (prot_de_mode == 'u'){
    uniprot_ids<-rownames(results_de)
    gene_symbols_all<-get_symbol_from_uniprot(uniprot_ids)
    gene_symbols<-gene_symbols_all$SYMBOL
   }else{
    gene_symbols<-rownames(results_de)
   }



  yvol = 'adj.P.Val'
     ## y = 'P.Value', 
    #results_de<-de_all_groups_proteins[[3]]
    pvol<-EnhancedVolcano(results_de, 
                    lab = gene_symbols,
                    x = 'logFC',
                    y = yvol,               
                    pCutoff = 0.05,
                    FCcutoff = 0.1, 
                    ylim=c(0,ylim), 
                    xlim=c(-xlim, xlim), 
                    title='', 
                    subtitle=ns
    )
    pvol
    fname_vol<-paste0(outdir_s_p,'/Volcano_', prefix_full,'_',VISIT_COMP,'_cluster_',cluster_id_num,'.jpeg')
    print(fname_vol)
    ggsave(fname_vol,pvol, width=6,height=8, dpi=300)
}

#ggsave(fname,pvol, width=6,height=8, dpi=300)

## Create a p-adjusted
## how many total proteins?
print(paste('Number of DE proteins in: ', tissue, VISIT_COMP))
de_proteins<-lapply(de_all_groups_proteins,function(x){(length(which(x$adj.P.Val<0.05)))})
print(unlist(de_proteins))

length(which(de_all_groups_proteins[[2]]$adj.P.Val<0.05))
dim(vsn_mat)[1]

#common_de<-intersect(all_sig_proteins,anova_results_oneway_significant)
#dir.create(outdir_s_p)
outdir_s_p_enrich<-paste0(outdir_s_p, '/enr_prot', tissue,'_', prot_de_mode, '/'); dir.create(outdir_s_p_enrich)
#write.csv(common_de, paste0(outdir_s_p, 'common_de.csv'))



#fit.cont_sig[common_de]
#



################### HEATMAPS  ############
# Heatmaps of the vs values OR the log2FC values per cluster 


#ARRANGE
#df_ord<-df[order(df$COHORT),]
pvol
order_by_metric<-'padj_reverse'
if (tissue=='CSF' & VISIT_COMP=='V08'){
  log2fol_T_overall=1
  
}
log2fol_T_overall=0
log2fol_T_hm=0
padj_T_overall_prot=0.05
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



################# ENRICHMENT - GSEA-GO #############

order_statistic<-'logFC'
#TODO: need to clean up remove gsea 

for (cluster_id_num in 1:4){

      results_de = de_all_groups_proteins[[cluster_id_num]]
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




      ################### run gsea with anova ######################
      gene_list1<-results_de[,order_statistic]
      names(gene_list1)<-rownames(results_de)
      gene_list_ord=gene_list1[order(-gene_list1)]

      
      
      
      pval_T_overall<-0.05
      gene_list_limma_significant=rownames(results_de)[results_de$adj.P.Val<padj_T_overall]
      gene_list_limma_significant_pval=rownames(results_de)[results_de$P.Value<pval_T_overall]
      gene_list_limma_significant_pval


      
      #gene_list_limma_significant

      run_anova=FALSE
      run_ORA=TRUE; 
      use_protein_pval=FALSE ## Proteins to use a sinput 
      use_pval=TRUE  ### WHAT TO PLOT## these do not work well with false so keep it and mark the number of sig in the legend


      pvalueCutoff_sig=0.05
   
      length(gene_list_ord)

      gene_list_ora=gene_list_limma_significant
      if (use_protein_pval){
        gene_list_ora=gene_list_limma_significant_pval
        
      }

      gene_list_ora
        gene_list_ord_sig<-gene_list_ord[names(gene_list_ord) %in% gene_list_ora]
          
          
        if (length(gene_list_ord_sig)>1){

        
         if (prot_de_mode == 'u'){
          uniprot_ids<-names(gene_list_ord_sig);
          uniprot_ids
          gene_symbols_all<-get_symbol_from_uniprot(uniprot_ids)
          gene_symbols_all
          gene_symbols<-gene_symbols_all$SYMBOL
          names(gene_list_ord_sig)<-gene_symbols
         }
          
    #  }else{
     ##     gene_symbols<-gene_list_ora
     # }
      ONT='BP'



      pvalueCutoff=0.05
      #outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT, '_', order_statistic)

      dir.create(paste0(outdir_s_p, '/enr_prot', tissue,'_', prot_de_mode,'/'), recursive=TRUE)
      outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, ONT,  '_', order_statistic, '_ora_', run_ORA, 'ppval_', use_protein_pval, 
                                    '_anova_', run_anova, 'pval_', use_pval, 'cl_', cluster_id_num )
      res_path<-paste0(outdir_s_p_enrich_file)

      outdir_s_p_enrich_file<-paste0(outdir_s_p_enrich, 'cl', cluster_id_num)
      outdir_s_p_enrich_file

      run_ORA=TRUE
      force_gse=TRUE
      # TODO: run per cluster 
      # Input:  gene_list_ord per cluster (coming from DE limma results )
      # 
              if (run_ORA){
                # writes 
                if (file.exists(paste0(outdir_s_p_enrich_file, '.Rds'))| !force_gse){

                      print('skip')
                }
              else{
#gene_list_ord_sig
#gene_list_ord_sig
                 if (length(gene_list_ord_sig)>0){
                  print('Running gse proteins')
                gse_protein_full = run_ora_gene_list(  gene_list_ord = gene_list_ord_sig, results_file_ora =outdir_s_p_enrich_file  )
                  saveRDS(gse_protein_full, paste0(outdir_s_p_enrich_file))
                 }

              }
              
        }
        
      #}outdir_s_p_enrich_file

        }

}




