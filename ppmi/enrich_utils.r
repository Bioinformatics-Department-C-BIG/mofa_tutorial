

# Cluster profiler enrichment analysis utils
calculate_log_fcs<-function( gse_all_cls, clust_names = c('BL', 'V06', 'V08'), metric='logfc'){
    # for a set of gene lists and their logFCs 
    # calculate average logFC
    #' @param gse_all_cls the gse_compare list of results
    #' @param clust_names
    #' 
    clust_names<-names(gse_all_cls)
    #gse_id=1
    log_fcs_all_clusts<-lapply(clust_names,function(gse_id ){
        
        # Calculate the logFCs
        gene_list_cluster_1<-gse_compare_cl@geneClusters[[gse_id]]
        gse_cluster_1<-gse_all_cls[[gse_id]] # select the results to obtain the DE pathways 

        colnames(gse_cluster_1)
       
        paths <- gse_cluster_1$Description # which are DE? 
        genes_in_path<-strsplit(gse_cluster_1$core_enrichment, '/')   # filter the significant? 

        
        # logFC per pathway 
        log_fcs<-lapply(genes_in_path,function(genes){
                return(mean(gene_list_cluster_1[genes], na.rm=TRUE))

        } )
        log_fcs_vec<-data.frame(logFC=unlist(log_fcs))

        #print(length(paths))
        log_fcs_vec$Description<-unlist(paths)
        log_fcs_vec<- cbind(log_fcs_vec, gse_cluster_1)

       # print(rownames(log_fcs_vec))
        return(log_fcs_vec)
      }
      )

    return(log_fcs_all_clusts)

  



}




   

    get_pathway_metrics_df<-function(log_fcs_all_clusts_list,clust_names,  metric='NES',padjust_cutoff=0.05 ){
    #'
    #' creates a merged result for all time points for a clusts 
    #' merges by pathway
    #' @param log_fcs_all_clusts_list holds all time ponts list of gse results 
    #' @param
    #' 
  #log_fcs_all_clusts_list = log_fcs_all_tps_list

      log_fcs_all_clusts_list_nes = lapply(log_fcs_all_clusts_list, function(x1){
        x2<-  x1[x1$p.adjust<padjust_cutoff,]

        x2[,c('Description', metric)]
      })

  
      # merge time points 
      merged_df<-Reduce(function(x, y) merge(x, y,  by='Description', all=TRUE), log_fcs_all_clusts_list_nes)
      #rownames(merged_df)<-merged_df$Description; merged_df$Description=NULL
      colnames(merged_df)<-c('Description' ,clust_names)

      return(merged_df)


}







#pvalueCutoff_mofa_factors=0.05
#factor=23

#fact # get factors for metric 
# get pathways for factors!! 

get_pathways_for_factors<-function(outdir, factors, pvalueCutoff_mofa_factors = 0.05, top_ps=10){
    #'Returns a set of significant pathways for factors
    #'@param oudtir
        gse_pathways_fs_list=lapply(factors, function(factor){
            results_file_mofa_f = paste0(outdir, '/enrichment/gsego_',factor,'_',pvalueCutoff_mofa_factors, '.csv')
            gse_factor_pathways<-read.csv(results_file_mofa_f)
            gse_factor_pathways<-gse_factor_pathways %>% dplyr::filter(p.adjust<0.05)
            gse_factor_pathways<-gse_factor_pathways %>% dplyr::arrange(p.adjust)
            print(head(gse_factor_pathways[,c('Description', 'p.adjust')]))
            return_fs_paths<-data.frame(Description=gse_factor_pathways$Description[1:top_ps])
            return_fs_paths$factor=factor
            return(return_fs_paths)
        })
        return(gse_pathways_fs_list)




}
