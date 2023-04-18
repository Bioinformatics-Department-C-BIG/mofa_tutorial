

source(paste0(script_dir, '/RNAseq enrichment.R'))
#### Prerequisites are mofa analysis 
### and rna seq enrichment 

list1=list()
nfactors<-MOFAobject@dimensions$K
mofa_enrich_rds<-paste0(outdir, '/enrichment/gse_results_mofa.Rds')

ONT='BP'

if (file.exists(mofa_enrich_rds)){
  list1<-loadRDS(mofa_enrich_rds)
  
}else{
  
      for (factor in 1:nfactors){
      
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
      
      saveRDS(list1,mofa_enrich_rds )

      
}


#gse_mofa=list1[[factor]]

## or LOAD GSE RESULTS HERE 

for (factor in 1:nfactors){
    
    
    results_file_mofa = paste0(outdir, '/enrichment/gsego_',factor,'_')
    gse_mofa=list1[[factor]]
    write.csv(as.data.frame(gse@result), paste0(results_file, '.csv'))
    
    ### to run mofa results
    run_enrichment_plots(gse=gse_mofa, results_file = results_file_mofa)
    
    
   
      
      
    }
    
    