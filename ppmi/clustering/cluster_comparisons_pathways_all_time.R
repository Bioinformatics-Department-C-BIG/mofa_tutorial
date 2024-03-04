

source(paste0(script_dir, 'ppmi/enrich_utils.R'))

  get_de_results_path<-function(deseq_params_all=deseq_params_all, VISIT=VISIT, formula_deseq_format=formula_deseq_format,  prefix=prefix, cluster_id=cluster_id ){
    #' returns the deseq results pathway for 
    #' @param VISIT
    #' @param formula_deseq_format
    #' @param prefix
    #' @param cluster_id
    #' 
            return(paste0(deseq_params_all,'/', VISIT, '/' ,formula_deseq_format, '/', prefix, 'de_cluster_', cluster_id , '.csv'))
        }

#Apply to current combinations 
# Clust comparisons for RNA/mirna

#' Enrichment settings for RNA/miRNA
#' This script takes the enrichment analysis results for all time points and concatenates them 
#' 
#' @value Produces a heatmap
#' Compares the de pathways 
order_by_metric='log2FoldChange'; order_by_metric_s='log2FC'
pvalueCutoff_sig=0.05



clusters_indices=c('1','2','3')
dir.create(deseq_params_all, '/all_time/enr/')
    ### Cluster compare by visit ### 
    # 1. load all gene lists again
    # deseq2ResDF<-read.csv(paste0(de_file), row.names=1 )
#formula_deseq_format
  if (!cell_corr_deseq){
    formula_deseq_format=''
  }

    times_sel<-c('BL', 'V06', 'V08')
    gse_compare_all_vis <-list()
     cluster_id='1'


    force_compare_time=FALSE
    for (cluster_id in clusters_indices){

      # cluster compare all visits together  
    #print(colnames(deseq2ResDF_time))
     # deseq_all_times<-vector("list", length = 3)

      deseq_all_times<-sapply(times_sel, function(VISIT){
        
        # deseq_params_all already loaded in cluster comparison script 
        de_results_path<-get_de_results_path(deseq_params_all, VISIT, formula_deseq_format,  prefix, cluster_id)

        deseq2ResDF_time<-read.csv(de_results_path, row.names=1) 
        
        gene_list1<-get_ordered_gene_list(deseq2ResDF_time,  order_by_metric, padj_T=1, log2fol_T=0 )
        names(gene_list1)<-gsub('\\..*', '',names(gene_list1))
        
        return(gene_list1)
        }
      )

      # Cluster compare for each cluster - compare the time points 
      dir.create(paste0(deseq_params_all, '/all_time/enr/', formula_deseq_format), recursive = TRUE)
      enrich_compare_path=paste0(deseq_params_all, '/all_time/enr/', formula_deseq_format, '/', prefix, enrich_params, cluster_id, 'time')
      if (!file.exists(paste0(enrich_compare_path, '.Rds' )) | force_compare_time){
     # if (TRUE){

      
            gse_compare_visit<-compareCluster(geneClusters = deseq_all_times, 
                                        fun = "gseGO", 
                                        OrgDb='org.Hs.eg.db', 
                                        ont=ONT, 
                                        keyType = 'ENSEMBL', 
                                        pvalueCutoff=1)  # allow them all in to get log2FC for all paths 


            plot_enrich_compare(gse_compare_visit,paste0(enrich_compare_path,cluster_id), N_EMAP = 60, N_DOT=7, N_DOT_U=10)

            saveRDS(gse_compare_visit,paste0(enrich_compare_path, '.Rds' ))

             gse_compare_all_vis[[cluster_id]]<-gse_compare_visit
            }else{
              gse_compare_all_vis[[cluster_id]]<-loadRDS(paste0(enrich_compare_path, '.Rds' ))
            }
  }

# TODO: load one with pvaluecutoff 1 and one with 0.05 





### LOAD ALL VISITS, ALL clusters and compare unique and intersection  ####

intersection_all_clusts=list()  # holds the intersection of pathways of all time points for each cluster  - length: nclusters 
union_all_clusts=list() # holds the union of pathways of all time points for each cluster  - length: nclusters 

for (cluster_id in clusters_indices){

      enrich_compare_path=paste0(deseq_params_all, '/all_time/enr/',formula_deseq_format,'/', prefix, enrich_params, cluster_id, 'time')

      # Print each visit 
      gse_compare_visit_res<-gse_compare_all_vis[[cluster_id]]@compareClusterResult
      # todo: FILTER by pvalue gse_compare_visit_res

      gse_compare_visit_res_t<-split(gse_compare_visit_res$Description,  gse_compare_visit_res$Cluster) # split by visit
      gse_compare_visit_res_t_sub<-gse_compare_visit_res_t 

      
       # Intersection of all time points
      intersection_clust<-Reduce( intersect,gse_compare_visit_res_t_sub)
      union_clust<-Reduce(union,gse_compare_visit_res_t_sub )
      length(unique(unlist(gse_compare_visit_res_t_sub)))

      length(union_clust)
      # save intersection and union for each cluster 
      intersection_all_clusts[[cluster_id]]<-intersection_clust
      union_all_clusts[[cluster_id]]<-union_clust

      # Venn diagram for each cluster and all time points 
      fname_venn=paste0(enrich_compare_path, 'venn_clust.png')
      create_venn(venn_list = gse_compare_visit_res_t_sub, fname_venn =fname_venn,main =paste0( 'Pathways by visit, cluster: ', cluster_id  ))


}

# [1] "Cluster"         "ID"              "Description"     "setSize"
# [5] "enrichmentScore" "NES"             "pvalue"          "p.adjust"
# [9] "qvalue"          "rank"            "leading_edge"    "core_enrichment"




# split for all clusters
cluster_id = '2'
sig_clust = 'V08'
padjust_cutoff = 0.05
padjust_cutoff = 0.001


#' @param sig_clust: select top paths from which time point? 
cluster_id = '2' ; sig_time = 'BL'; 

sig_time = 'V08'

get_top_per_clust<-function(gse_compare_all_vis,top_paths=top_paths,sig_time=sig_time ){

  # return the top significant pathways by cluster
  # sig_time: for a specific time point
  top_clust<-lapply(gse_compare_all_vis, function(gse_compare_cl){
    gse_all_cls<-split(gse_compare_cl@compareClusterResult,  gse_compare_cl@compareClusterResult$Cluster) # split by visit
    gene_list_cluster_1<-gse_compare_cl@geneClusters[[sig_time]]
    gse_cluster_1<-gse_all_cls[[sig_time]]
    
    gse_cluster_1<-gse_cluster_1 %>% filter(p.adjust<0.05)

    paths1_sig_top<-gse_cluster_1$Description[order(gse_cluster_1$p.adjust)]


    if (top_paths){
      paths1_sig_top<-paths1_sig_top[1:top_paths]
    }
    #paths1_sig<-gse_cluster_1$Description[gse_cluster_1$p.adjust<padjust_cutoff] # also cut with padjusted? 
    paths1_sig_top<-na.omit(paths1_sig_top) # only plot top 
    return(paths1_sig_top)


  })
  return(top_clust)
 
}



cluster_id = '2' ; sig_time = 'V08';
top_paths = 60

top_paths_all_factors<-concatenate_top_pathways_factors(fact, pvalueCutoff = 0.05, top_p = top_paths)
dim(top_paths_all_factors)
top_paths_all_factors$Description
#'  @param metric

metric='logFC';
gse_compare_cl=gse_compare_all_vis[[1]]
#' decide on pathways to keep 
top_sig_all_clusts
all_sig_all_clusts<-unlist(get_top_per_clust(gse_compare_all_vis, top_paths = FALSE,sig_time=sig_time))
top_sig_all_clusts<-unlist(top_sig_all_clusts)
metric='NES';
metric = 'logFC'

# holds all timepoints all, clusters 
#
padjust_hm<-1
log_fcs_all_tps_all_clusts<-lapply(gse_compare_all_vis, function(gse_compare_cl){
  #' @param 
  #' @param
  #' gse_compare_cl<-[[cluster_id]]
  #' @return merged_df_lfc a dataframe with logFC/NES for all timepoints (n paths x n timepoints)
  gse_all_cls<-split(gse_compare_cl@compareClusterResult,  gse_compare_cl@compareClusterResult$Cluster) # split by visit

  log_fcs_all_tps_list<-calculate_log_fcs(gse_all_cls ) # calculates the average logFC per pathway
  

  merged_df_lfc<-get_pathway_metrics_df(log_fcs_all_tps_list,metric=metric, clust_names=names(gse_all_cls) , padjust_cutoff = padjust_hm)

  return(merged_df_lfc)
})


#which_n<-log_fcs_all_tps_list[[3]]$Description =='natural killer cell mediated cytotoxicity'
#log_fcs_all_tps_list[[3]][which_n,]

log_fcs_all_tps_all_clusts<-lapply(1:length(log_fcs_all_tps_all_clusts), function(i){
  x = log_fcs_all_tps_all_clusts[[i]]
  list_names<-names(log_fcs_all_tps_all_clusts)
  colnames(x)<-paste0(colnames(x), '_', list_names[[i]]) 
  colnames(x)[1]<-'Description'
  return(x)
}
)



# merge all dfs for all clusters
merged_df_all_tps_all_clusts <-Reduce(function(x, y) merge(x, y,  by='Description', all=TRUE), log_fcs_all_tps_all_clusts)


# get a merged df per cluster 



length(all_sig_all_clusts)

logFC_merged_df2<-merged_df_all_tps_all_clusts
use_mofa_paths = FALSE

if (use_mofa_paths){
  top_paths = 30

  top_paths_all_factors<-concatenate_top_pathways_factors(fact, pvalueCutoff = 0.05, top_p = top_paths)
    dim(top_paths_all_factors)

  top_paths_all_factors<-top_paths_all_factors[top_paths_all_factors$Description %in% all_sig_all_clusts,]
   dim(top_paths_all_factors)
  selected_paths<-top_paths_all_factors$Description


}else{
  top_paths=15
  top_sig_all_clusts<-unlist(get_top_per_clust(gse_compare_all_vis,  top_paths = top_paths,sig_time=sig_time))

  selected_paths<-top_sig_all_clusts
}

length(top_sig_all_clusts)
# filter top 
  logFC_merged_df2<-logFC_merged_df2[logFC_merged_df2$Description %in% selected_paths,]

dim(logFC_merged_df2)
logFC_merged_df2[is.na(logFC_merged_df2)]<-0

rownames(logFC_merged_df2)<-logFC_merged_df2$Description
logFC_merged_df2$Description<-NULL
graphics.off()

enrich_compare_path_all_clusts=paste0(deseq_params_all, '/all_time/enr/', formula_deseq_format, '/', prefix, enrich_params)
enrich_compare_path_all_clusts

fname= paste0(enrich_compare_path_all_clusts,'clall_' , metric,'mofa_',as.numeric(use_mofa_paths),'_', top_paths, '_heatmap.png')
fname
row_an<-top_paths_all_factors[match( rownames(logFC_merged_df2), top_paths_all_factors$Description  ), ] 
row_an2<-as.factor(row_an$factor); names(row_an2)<-row_an$Description

row_ha = rowAnnotation( factor=row_an2)

png(fname, width=20*100, height=20*100, res=200)
ch<-ComplexHeatmap::pheatmap(as.matrix(logFC_merged_df2), 
    show_rownames=TRUE, 
    column_split = rep(1:3, each=3), 
    right_annotation = row_ha, 
    heatmap_legend_param  = list(direction = "horizontal"))

    
    draw(ch, heatmap_legend_side="bottom", padding = unit(c(2, 2, 2, 70), "mm")
    )

  #  title = paste(metric, 'p-adj:', padjust_hm))
dev.off()

logFC_merged_df2[10, ]





fname_venn=paste0(enrich_compare_path, 'time_all_intersection.png')
create_venn(venn_list = intersection_all_clusts, fname_venn =fname_venn,main =paste0( 'Pathways by cluster' ))

fname_venn=paste0(enrich_compare_path, 'time_all_union.png')
create_venn(venn_list = union_all_clusts, fname_venn =fname_venn,main =paste0( 'Pathways by cluster - union' ))



Reduce(intersect,intersection_all_clusts)
Reduce(intersect,union_all_clusts)


# Choose union or intersection to compare
unique_partitions<-VennDiagram::get.venn.partitions(union_all_clusts)

unique_partitions



library(tidyverse)

lists <- intersection_all_clusts
lists <- union_all_clusts

unique_partitions<-data.frame(data = names(lists), number = matrix(lists)) %>%
  tidyr::unnest(cols = c(number)) %>%
  dplyr::distinct(number, .keep_all = TRUE) %>% as.data.frame()

sapply(c('1','2','3'), function(x){
    unique_partitions_c<-unique_partitions[unique_partitions$data==x,]
    write.csv(unique_partitions_c,paste0(enrich_compare_path, '_union_unique',x,'.csv' ))

})


lists <- intersection_all_clusts

unique_partitions_inter<-data.frame(data = names(lists), number = matrix(lists)) %>%
  tidyr::unnest(cols = c(number)) %>%
  dplyr::distinct(number, .keep_all = TRUE) %>% as.data.frame()

sapply(c('1','2','3'), function(x){
    unique_partitions_inter_c<-unique_partitions_inter[unique_partitions_inter$data==x,]
    unique_partitions_union_c<-unique_partitions[unique_partitions$data==x,]

    write.csv(unique_partitions_inter_c,paste0(enrich_compare_path, '_intersection_unique',x,'.csv' ))
    write.csv(unique_partitions_union_c,paste0(enrich_compare_path, '_union_unique',x,'.csv' ))

    unique_partitions_i_un_c<-intersect(unique_partitions_inter_c, unique_partitions_union_c )

    write.csv(unique_partitions_union_c,paste0(enrich_compare_path, '_union_inter_unique',x,'.csv' ))

    



})



#' @gse_compare_visit_res_t: holds all times for each cluster
unique_timepoint<-list()
dir.create(paste0(deseq_params_all, '/all_time/enr/'))
sel_t<-'V08'
# check Swhat is in V08 only AND in the unique of the union of all times
for (cluster_id in clusters_indices){

      enrich_compare_path=paste0(deseq_params_all, '/all_time/enr/', prefix, enrich_params, cluster_id, 'time')


      gse_compare_visit_res<-gse_compare_all_vis[[cluster_id]]@compareClusterResult
      # split by cluster 
      names(gse_compare_visit_res)
      gse_compare_visit_res_t<-split(gse_compare_visit_res,  gse_compare_visit_res$Cluster) 

    unique_partitions_union_c<-unique_partitions[unique_partitions$data==cluster_id,]

    unique_timepoint[[cluster_id]]<-gse_compare_visit_res_t[[sel_t]][gse_compare_visit_res_t[[sel_t]]$Description %in% unique_partitions_union_c$number, ]
    gse_compare_visit_res_t[[sel_t]]

    print(dim(unique_timepoint[[cluster_id]]))

    # 
    write.csv(unique_timepoint[[cluster_id]],paste0(enrich_compare_path, 'unique',sel_t,'_',cluster_id,'.csv' ))


}
gse_compare_visit_res_t


enrich_compare_path




## is the intersection also in the unique of the union? 
### PLOTTING 

## gse
## Plot the enrichment score for each cluster over time


gse_clust_pathway=list()







for (cluster_id in clusters_indices){
    print(cluster_id)

      enrich_compare_path=paste0(deseq_params_all, '/enr/', prefix, enrich_params, cluster_id, 'time')


      gse_compare_visit_res<-gse_compare_all_vis[[cluster_id]]@compareClusterResult
      # split by time
      gse_compare_visit_res_t<-split(gse_compare_visit_res,  gse_compare_visit_res$Cluster) 

      length(gse_compare_visit_res_t)

      ## apply for each visit, take the enrichment score and plot it 
      cols<-c('Description', 'p.adjust')
    gse_clust_pathway[[cluster_id]]<-lapply(gse_compare_visit_res_t, function(gse_clust){
      return(gse_clust[grep('MHC', gse_clust$Description, ), cols])
    })



}


gse_clust[grep('lipo', gse_clust$Description ),'Description']


gse_clust_pathway[[1]]




clust_id=3
dim(gse_clust_pathway[[clust_id]][[2]])
np<-dim(gse_clust_pathway[[clust_id]][[1]])[2];np
new<-do.call(rbind,gse_clust_pathway[[clust_id]][c(1,2,3)])

new$VISIT<-c(rep('V06',np), rep('V08', np))
new$VISIT

ggplot(new,aes(x=VISIT, y=-log10(p.adjust), group=Description))+
    geom_line(aes(color=Description))










































































