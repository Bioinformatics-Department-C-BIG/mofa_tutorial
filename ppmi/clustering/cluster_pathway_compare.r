
## This script aims to find unique pathways in 1,2,3
clust_comparison_indices<-c('1',  '1_2')
to_remove_clust<-c()

clust_comparison_indices<-c('1','1_2',  '1_3')
to_remove_clust<-c('2','3')
to_remove_clust<-c('2','3')
to_remove_clust<-c()

# TODO: 3,2
clust_id = '1'


#clust_comparison_indices<-c('3','3_1',  '3_2')
#to_remove_clust<-c('1','2')
#clust_id = '3'



gse_compare_df<-as.data.frame(gse_compare)
gse_compare_sig<-gse_compare_df[gse_compare_df$p.adjust<0.05,]

gse_compare_sig_by_clust<-split(gse_compare_sig$Description, gse_compare_sig$Cluster)
gse_compare_sig_by_clust_pv<-split(gse_compare_sig, gse_compare_sig$Cluster)
names(gse_compare_df_sig_by_clust)
# gse_compare_sig: holds the pathway descriptions for the DE in cluster 1

# intersection of the ones DE with HC AND other clusters 
gse_compare_sig_inter<- Reduce( intersect,gse_compare_sig_by_clust[clust_comparison_indices])
gse_compare_sig_inter


# remove the ones DE with control and other clusters 
if (length(to_remove_clust)>0){
gse_compare_sig_inter_unique<-gse_compare_sig_inter[!(gse_compare_sig_inter %in% unlist(gse_compare_sig_by_clust[to_remove_clust]))]

}else{
    gse_compare_sig_inter_unique = gse_compare_sig_inter
}



# Print the 3 lists 

clust1<-gse_compare_sig_by_clust_pv[[clust_id]]

gse_compare_sig_inter_pvals<-clust1[clust1$Description %in% gse_compare_sig_inter_unique,c('Description','p.adjust') ]
gse_compare_sig_inter_pvals

paste0(enrich_compare_path, paste0(clust_comparison_indices, collapse='-'), '.csv')
write.csv(gse_compare_sig_inter_pvals,paste0(enrich_compare_path, paste0(clust_comparison_indices, collapse='-'), '.csv') )


print(gse_compare_sig_inter_pvals)















