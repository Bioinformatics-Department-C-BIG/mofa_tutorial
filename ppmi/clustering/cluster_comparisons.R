


# WARNING DO NOT CHANGE THIS, MOVE TO A FUNCTION!! ## 
### FIRST LOAD required files  ####
source(paste0(script_dir, '/ppmi/utils.R'))

process_mirnas = TRUE; # reload mirs !!  # DO NOT CHANGE THIS!! 
source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se_mirs=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
se_mirs_norm=load_se_all_visits(input_file = input_file_mirs, combined=combined_bl_log); 
se_mirs
head(assay(se_mirs))
head(assay(se_mirs_norm))

process_mirnas=FALSE
source(paste0(script_dir, '/ppmi/config.R'))
#source(paste0(script_dir, '/ppmi/deseq2_vst_preprocessing_mirnas_all_visits2.R'))
se_rnas=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 


## 1. get Summarized Experiment with metrics from all time points 
## 2. Run deseq 
## 3. enrichment 

MOFAobject_clusts=MOFAobject_sel # take it from the clusterig of the last visit only 


## SETUP the script parameters here ####
# 1. select visit, 2. process mirs 
# TODO: make function to load for rnas and mirnas separately
# edit this one 
process_mirnas= TRUE
if (process_mirnas){
  se_sel = se_mirs
  prefix='mirnas_'
}else{
  se_sel = se_rnas
  prefix='rnas_'
}

view=ifelse(process_mirnas, 'miRNA', 'RNA');view


se_clusters<-filter_se(se_sel, VISIT=VISIT_COMP, sel_coh = sel_coh, sel_sub_coh = sel_ps)
se_clusters
### Decide on the parameters settings 
# Set the outdirectory 

y_clust='moca'
clust_name=paste0(y_clust, '_clust')

## Outputs 
# 1. DE files 
# 2. Venns 
# 3. Volcano plot
# 4. Enrichment analysis 


#MOFAobject_clusts<-MOFAobjectPD
deseq_all_groups <- vector("list", length = 3);
se_filt_all<- vector("list", length = 3);

# Correct for blood cell proportions of neutrophils and lymphocytes 
cell_corr<-TRUE

# Obtain clustering from mofa
se_clusters$kmeans_grouping<- groups_from_mofa_factors(se_clusters$PATNO, MOFAobject_clusts, y_clust );
se_clusters$kmeans_grouping=as.numeric(se_clusters$kmeans_grouping)
nclusts = length(table(se_clusters$kmeans_grouping));nclusts
#cluster_params_dir<-paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/')
cluster_params<-paste0(clust_name ,'/', k_centers_m,'/r',as.numeric(rescale_option),'/g', as.numeric(sel_group_cors) )
cluster_params_dir<-paste0(outdir,'/clustering/',cluster_params );
cluster_params_dir

cd <- colData(se_clusters)
colData(se_clusters)[cd$INEXPAGE%in%'INEXHC','kmeans_grouping']<-'HC'
se_clusters$kmeans_grouping<-as.factor(se_clusters$kmeans_grouping)




# todo move to config
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
formula_deseq
# TODO: include corrected for cell counts  and uncorrected formula
#if (process_mirnas){
  # TODO: try site? and lymphocytes too? 
#  formula_deseq = '~AGE_SCALED+SEX+kmeans_grouping'

#  cell_corr = FALSE # ensure they are always in the same folder c0
# apply this both for mirs and rnas 
if (cell_corr) {
  formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+Neutrophil.Score+Plate+kmeans_grouping'

}else{
  formula_deseq = '~AGE_SCALED+SEX+Plate+Usable_Bases_SCALE+Plate+kmeans_grouping'


}
formula_deseq
cluster_id=1

#param <- SnowParam(workers = 6, type = "MPI")


deseq_all <- vector("list", length = 0) # holds the significant gene/mirs ids only for each cluster
deseq_all_names<-vector("list", length = 0)
gene_lists<-vector("list", length = 0)
gse_all<-vector("list", length = 0)

### se_filt_all: list to hold the se
### deseq_all_groups: list to hold the deseq results 
### deseq_significant_all_groups: list to hold significant 

deseq_params_all<-paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr))
deseq_params<-paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr),  '/',VISIT_COMP)
dir.create(paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr)))
dir.create(deseq_params)
fname_venn=paste0(deseq_params,'/', prefix , 'min_',min.count,'venn.png');fname_venn

deseq_params

# TODO: ensure thatthis is also run for all subjects together
for (cluster_id in 1:3){

  ### 1. for each cluster, create se filt with controls, 
  ### 2. run deseq 
  ### 3. get significant per cluster 
  print(paste('cluster:',cluster_id))

  #de_file
  se_filt_all[[cluster_id]]<-se_clusters[,se_clusters$kmeans_grouping %in% c(cluster_id,'HC')]
}
## filter for grouping here - it is not in the prerpocessing by default
#se_filt<-se_filt[,!(is.na(se_filt$kmeans_grouping))] 

formula_deseq
se_filt_all

clusters_indices= list( '1','2','3', c(1,2,3))
clusters_indices
for (cluster_id in clusters_indices){

  ### 1. for each cluster, create se filt with controls, 
  ### 2. run deseq 
  ### 3. get significant per cluster 
  #print(paste('cluster:',cluster_id))
  cluster_id_index=paste0(cluster_id, collapse='_')
  de_file<-paste0(deseq_params, '/', prefix, 'de_cluster_', paste0(cluster_id, collapse='_') , '.csv')
  #de_file
  se_filt_all[[cluster_id_index]]<-se_clusters[,se_clusters$kmeans_grouping %in% c(cluster_id,'HC')]


   if (file.exists(de_file)){
  #if (FALSE)
    # if de file exists load it - unfiltered de results file
    deseq2ResDF<-read.csv(paste0(de_file), row.names=1 )

  }else{
    # else run the deseq with the design formula specified 
        deseq2ResDF = deseq_by_group(se_filt_all[[cluster_id_index]], formula_deseq, min.count=min.count)

        deseq_all_groups[[cluster_id_index]]<-deseq2ResDF
        if (!process_mirnas){
          # get symbols for RNA only 
          deseq2ResDF$GENE_SYMBOL<-get_symbols_vector(gsub('\\..*', '',rownames(deseq2ResDF))) 
        }
        write.csv(deseq2ResDF, de_file, row.names=TRUE)
  }
  print(cluster_id_index)
  deseq_all_groups[[cluster_id_index]]<-deseq2ResDF
  deseq_all[[cluster_id_index]]<-deseq2ResDF[deseq2ResDF$mofa_sign %in% 'Significant',] # holds the significant only
} 
# Save and load # Rrename ens id.*


deseq_all_names <- lapply(deseq_all, function(x){return(  gsub('\\..*', '',rownames(x))   )  })
names(deseq_all_names) <-names(deseq_all)




#### 1. Venn from significant 
length(deseq_all_names)
create_venn(venn_list = deseq_all_names, fname_venn =fname_venn,
            main =paste0( ' DE molecules for each molecular cluster' ))

graphics.off()
 # TODO: 

########### 



# TODO: venn before and after correction 

for (cluster_id in clusters_indices){

      # Take the original data and deseq results to plot a volcano 
      print(paste('cluster:',cluster_id))
      cluster_id_index=paste0(cluster_id, collapse='_')
      se_filt=se_filt_all[[cluster_id_index]]

      deseq2ResDF=deseq_all_groups[[cluster_id_index]]
      deseq2ResDF$log2FoldChange
      #deseq2ResDF$GENE_SYMBOL
      deseq2ResDF$padj

      pvol<-plotVolcano(  deseq2ResDF, se_filt, title=paste0('Cluster ', cluster_id_index), xlim=c(-1.1,1.1),
      lab=deseq2ResDF$GENE_SYMBOL)
      pvol
      fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT_COMP,'.jpeg')
      fname<-paste0(deseq_params, '/Volcano_', prefix, VISIT_COMP,'_cluster_', cluster_id_index, '.jpeg')

      pvol
      ggsave(fname,pvol, width=9,height=12, dpi=300)
}


#### 2. Venn from significant in top of factor

## 1. get top of factor / or all higly variable genes input into MOFA
## 2. intersect with DE 
sel_factor=4
top10<- gsub('\\..*', '',select_top_bottom_perc(MOFAobject=MOFAobject, view=view, factors = sel_factor, top_fr = 0.1))
top10<- gsub('\\..*', '',select_top_bottom_perc(MOFAobject=MOFAobject, view=view, factors = sel_factor, top_fr = 1))
length(top10)
# intersect with the top factors 
deseq_all_top<-lapply(deseq_all_names, function(x) intersect(x,top10 ) )
deseq_all_top


fname_venn=paste0(deseq_params, prefix, 'venn_de_per_group_deseq', 'top_f', sel_factor ,  '.png')
create_venn(venn_list = deseq_all_top, fname_venn =fname_venn,main =paste0( ' DE molecules for each molecular cluster AND highly variable' ))





deseq2ResDF$log2FoldChange
order_by_metric='log2pval';order_by_metric_s='log2p'
order_by_metric='log2FoldChange'; order_by_metric_s='log2FC'

ONT='BP'
pvalueCutoff_sig=0.05
enrich_params<-paste0(ONT, order_by_metric_s)
dir.create(paste0(deseq_params, '/enr/'))

enrich_compare_path=paste0(deseq_params, '/enr/', prefix, enrich_params, 'comp')

results_file_cluster

for (cluster_id in clusters_indices){
  # run enrichment with the log2foldchange metric
  print(cluster_id) 

  cluster_id_index=paste0(cluster_id, collapse='_')
  deseq2ResDF = deseq_all_groups[[cluster_id_index]]
  gene_list1<-get_ordered_gene_list(deseq2ResDF,  order_by_metric, padj_T=1, log2fol_T=0 )
  names(gene_list1)<-gsub('\\..*', '',names(gene_list1))

  gene_lists[[cluster_id_index]]<-gene_list1 # load the gene lists separately for cluster compare 

  results_file_cluster=paste0(deseq_params, '/enr/', prefix, enrich_params, 'cl', cluster_id_index)
  gse_file<-paste0(results_file_cluster, '.Rds')
  if (!file.exists(gse_file)){
      gse1<-run_enrich_per_cluster(deseq2ResDF, results_file_cluster,N_DOT=20, 
      N_EMAP=30 , N_NET=10)
      saveRDS(gse1,gse_file )

  }
}


# Run cluster compare by cluster - it does not need the separate files only the gene lists 
clust_pair<-c(1,2,3)
clust_pair_s=paste0(clust_pair, collapse='_')
clust_pair_s
geneClusters = list(G1=gene_lists[[clust_pair[1]]],
                            G2=gene_lists[[clust_pair[2]]],  # nolint
                            G3=gene_lists[[clust_pair[3]]])
geneClusters=gene_lists
gene_lists
gse_compare<-compareCluster(geneClusters = geneClusters , 
                            fun = "gseGO", 
                            OrgDb='org.Hs.eg.db', 
                          ont=ONT, 
                          keyType = 'ENSEMBL') 


### RUN SCRIPT compare
names(geneClusters)                  
plot_enrich_compare(gse_compare,paste0(enrich_compare_path,clust_pair_s), N_EMAP = 60, N_DOT=8)




### Cluster compare by visit ### 
# 1. 
#    deseq2ResDF<-read.csv(paste0(de_file), row.names=1 )

cluster_id=3
deseq_all_times<-vector("list", length = 3)
deseq_all_times


deseq_all_times<-sapply( c('V04', 'V06', 'V08'), function(VISIT){
  deseq2ResDF_time<-read.csv(paste0(deseq_params_all,'/', VISIT, '/' , prefix, 'de_cluster_', cluster_id , '.csv'), row.names=1) 
  
  gene_list1<-get_ordered_gene_list(deseq2ResDF_time,  order_by_metric, padj_T=1, log2fol_T=0 )
  names(gene_list1)<-gsub('\\..*', '',names(gene_list1))


  return(gene_list1)
  }
)

dir.create(paste0(deseq_params_all, '/enr/'))
enrich_compare_path=paste0(deseq_params_all, '/enr/', prefix, enrich_params, cluster_id, 'time')


## Compare the three clusters for one visit 
enrich_compare_path
gse_compare<-compareCluster(geneClusters = list(T1=deseq_all_times[[1]],T2=deseq_all_times[[2]],T3=deseq_all_times[[3]] ), 
                            fun = "gseGO", 
                            OrgDb='org.Hs.eg.db', 
                            ont=ONT, 
                            keyType = 'ENSEMBL') 


plot_enrich_compare(gse_compare,paste0(enrich_compare_path,clust_pair_s), N_EMAP = 80)



























































































































