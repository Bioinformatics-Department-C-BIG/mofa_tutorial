





fact=get_factors_for_metric(y_clust)


top_paths_all_factors<-concatenate_top_pathways_factors(fact, pvalueCutoff = 0.05, top_p = top_paths, prefix='mirnas_')

top_paths_all_factors



# TODO: collect pathways by cluster for all time points 

results_file_cluster
enrich_params_mirs
deseq_params
get_enrich_result<-function(vis='V08', cluster_id_name = 2){

    results_file_cluster=paste0(deseq_params,'/../',vis,  '/enr_',prefix, '/', prefix, enrich_params_mirs, 'cl', cluster_id_name)
    results_file_cluster = paste0(results_file_cluster,'_', pvalueCutoff, '.csv' )
    return(results_file_cluster)
}

vis='BL'

all_times_paths<-list()
all_times_paths <- lapply(times_sel, function(vis){
    fname<-get_enrich_result(vis=vis)
    # IF THE FILE is there and there are more than 1 rows 
    if (file.exists(fname) && length(count.fields(fname, sep = ","))>1 ){
         return(read.csv2(fname,sep=','))
         

    }

}


)

is.null(all_times_paths)
colnames(all_times_paths[[3]])
all_times_paths[[3]]$Subcategory %in% top_paths_all_factors$Description

# TODO: get the log FCs of the mirnas - need to extract from the LOG FC input  files?
# supply the mir logFC  list? 
#merged_df_all_tps_all_clusts <-Reduce(function(x, y) merge(x, y,  by='Description', all=TRUE), log_fcs_all_tps_all_clusts)



fname









