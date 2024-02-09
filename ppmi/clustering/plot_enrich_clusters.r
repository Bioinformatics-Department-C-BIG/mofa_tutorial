



names(gse_all_clusters)


## gse
## Plot the enrichment score for each cluster over time
for (cluster_id in clusters_indices){
    print(paste('cluster', cluster_id))
    gse_clust<-gse_all_clusters[[cluster_id]]@result
    gse_clust$enrichmentScore
    print(gse_clust[grep('MHC', gse_clust$Description, ), c('Description','NES', 'enrichmentScore', 'p.adjust')])
}






