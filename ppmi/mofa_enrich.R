

list1=list()
num_factors<-MOFAobject@dimensions$K
for (factor in 1:num_factors){
      ONT='BP'

      f1<-get_weights(MOFAobject, view='RNA', factor=factor)
      gene_list<-f1['RNA'][[1]] 
      hist(gene_list)
      
      order_ind<-order(-gene_list)
      gene_list_ord<-gene_list[order_ind,]
      names(gene_list_ord)<-rownames(gene_list)[order_ind]
      names(gene_list_ord)

      gse_mofa <- clusterProfiler::gseGO(gene_list_ord, 
                                  ont=ONT, 
                                  keyType = 'ENSEMBL', 
                                  OrgDb = 'org.Hs.eg.db', 
                                  pvalueCutoff  = 0.05)
      
    list1[[factor]]<-gse_mofa

}
saveRDS(list1, paste0(outdir, 'gse_results_mofa.csv'))
saveRDS(list1, paste0(outdir, 'gse_results_mofa.csv'))




res1<-list1[[factor]]@result

View(gse@result)
run_mofa=TRUE

