

#Apply to current combinations 
# 

clusters_indices=c('1','2','3')


    ### Cluster compare by visit ### 
    # 1. load all gene lists again
    # deseq2ResDF<-read.csv(paste0(de_file), row.names=1 )

  if (!cell_corr_deseq){
    formula_deseq_format=''
  }

    times_sel<-c('BL', 'V06', 'V08')
    gse_compare_all_vis <-list()

    for (cluster_id in clusters_indices){

      # cluster compare all visits together  
    
      deseq_all_times<-vector("list", length = 3)

      deseq_all_times<-sapply(times_sel, function(VISIT){
        deseq2ResDF_time<-read.csv(paste0(deseq_params_all,'/', VISIT, '/' ,formula_deseq_format, '/', prefix, 'de_cluster_', cluster_id , '.csv'), row.names=1) 
        
        gene_list1<-get_ordered_gene_list(deseq2ResDF_time,  order_by_metric, padj_T=1, log2fol_T=0 )
        names(gene_list1)<-gsub('\\..*', '',names(gene_list1))
        
        return(gene_list1)
        }
      )





      # Cluster compare for each cluster - compare the time points 
      dir.create(paste0(deseq_params_all, '/enr/', formula_deseq_format))
      enrich_compare_path=paste0(deseq_params_all, '/enr/', formula_deseq_format, '/', prefix, enrich_params, cluster_id, 'time')
      if (!file.exists(paste0(enrich_compare_path, '.Rds' ))){
     # if (TRUE){

      
            gse_compare_visit<-compareCluster(geneClusters = deseq_all_times, 
                                        fun = "gseGO", 
                                        OrgDb='org.Hs.eg.db', 
                                        ont=ONT, 
                                        keyType = 'ENSEMBL') 


            plot_enrich_compare(gse_compare_visit,paste0(enrich_compare_path,cluster_id), N_EMAP = 80, N_DOT=5)

            saveRDS(gse_compare_visit,paste0(enrich_compare_path, '.Rds' ))

             gse_compare_all_vis[[cluster_id]]<-gse_compare_visit
            }else{
              gse_compare_all_vis[[cluster_id]]<-loadRDS(paste0(enrich_compare_path, '.Rds' ))
            }
  }







### LOAD ALL VISITS, ALL clusters and compare unique and intersection  ####

intersection_all_clusts=list()  # holds the intersection of pathways of all time points for each cluster  - length: nclusters 
union_all_clusts=list() # holds the union of pathways of all time points for each cluster  - length: nclusters 



gse_compare_all_vis
for (cluster_id in clusters_indices){

      enrich_compare_path=paste0(deseq_params_all, '/enr/', prefix, enrich_params, cluster_id, 'time')

      # Print each visit 
      gse_compare_visit_res<-gse_compare_all_vis[[cluster_id]]@compareClusterResult
      gse_compare_visit_res_t<-split(gse_compare_visit_res$Description,  gse_compare_visit_res$Cluster) # split by visit
      gse_compare_visit_res_t_sub<-gse_compare_visit_res_t 


    # Get
    


      
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

formula_deseq_format
# TODO: get the 
colnames(gse_compare@compareClusterResult)
# [1] "Cluster"         "ID"              "Description"     "setSize"
# [5] "enrichmentScore" "NES"             "pvalue"          "p.adjust"
# [9] "qvalue"          "rank"            "leading_edge"    "core_enrichment"




# split for all clusters


cluster_id = '2'

sig_clust = 'V08'
padjust_cutoff = 0.05


cluster_id = '1'
sig_clust = 'V08'
padjust_cutoff = 0.001

gse_compare_cl<-gse_compare_all_vis[[cluster_id]]
gse_all_cls<-split(gse_compare_cl@compareClusterResult,  gse_compare_cl@compareClusterResult$Cluster) # split by visit
gene_list_cluster_1<-gse_compare_cl@geneClusters[[sig_clust]]


gse_cluster_1<-gse_all_cls[[sig_clust]]
paths1_sig<-gse_cluster_1$Description[gse_cluster_1$p.adjust<padjust_cutoff]
length(paths1_sig)
gse_compare_cl

clust_nams<-names(gse_compare_cl@geneClusters)
names(gse_compare_cl)
# by cluster 
#c('1', '2', '3')
log_fcs_all_clusts<-lapply(c(1,2,3),function(clust_id){

        gene_list_cluster_1<-gse_compare_cl@geneClusters[[clust_id]]
        gse_cluster_1<-gse_all_cls[[clust_id]]
        paths<-gse_cluster_1$Description
        # filter the significant? 
        genes_in_path<-strsplit(gse_cluster_1$core_enrichment, '/')

        # logfc per pathway 
        log_fcs<-lapply(genes_in_path,function(genes){
                return(mean(gene_list_cluster_1[genes], na.rm=TRUE))

        } )
        log_fcs_vec<-as.data.frame(unlist(log_fcs))
        
        print(length(paths))
        log_fcs_vec$Description<-unlist(paths)
       # print(rownames(log_fcs_vec))
        return(log_fcs_vec)
}
)


merged_df<-Reduce(function(x, y) merge(x, y,  by='Description', all=TRUE), log_fcs_all_clusts)
rownames(merged_df)<-merged_df$Description; merged_df$Description=NULL
colnames(merged_df)
colnames(merged_df)<-clust_nams
merged_df2<-merged_df
merged_df2[merged_df2==NA]<-0
#merged_df2<-na.omit(merged_df)
#merged_df2<-merged_df
rownames(merged_df2) %in% paths1_sig
merged_df2<-merged_df2[rownames(merged_df2) %in% paths1_sig,]

fname= paste0(enrich_compare_path,'cl_',cluster_id, p.adjust_cutoff, '_heatmap.png')
jpeg(fname, width=10*100, height=15*100, res=150)

pheatmap(as.matrix(merged_df2))
dev.off()


paths
unlist(log_fcs)

gse_compare_visit_res_t_sub




length(intersection_all_clusts)

intersection_all_clusts

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
sel_t<-'V08'
# check Swhat is in V08 only AND in the unique of the union of all times
for (cluster_id in clusters_indices){

      enrich_compare_path=paste0(deseq_params_all, '/enr/', prefix, enrich_params, cluster_id, 'time')


      gse_compare_visit_res<-gse_compare_all_vis[[cluster_id]]@compareClusterResult
      # split by cluster 
      names(gse_compare_visit_res)
      gse_compare_visit_res_t<-split(gse_compare_visit_res,  gse_compare_visit_res$Cluster) 

    unique_partitions_union_c<-unique_partitions[unique_partitions$data==cluster_id,]

    unique_timepoint[[cluster_id]]<-gse_compare_visit_res_t[[sel_t]][gse_compare_visit_res_t[[sel_t]]$Description %in% unique_partitions_union_c$number, ]
    gse_compare_visit_res_t[[sel_t]]

    print(dim(unique_timepoint[[cluster_id]]))

    write.csv(unique_timepoint[[cluster_id]],paste0(enrich_compare_path, '_unique_per_timepoint',sel_t,'_',cluster_id,'.csv' ))


}
gse_compare_visit_res_t







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








































































