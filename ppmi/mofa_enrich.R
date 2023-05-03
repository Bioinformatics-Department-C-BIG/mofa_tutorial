

source(paste0(script_dir, '/RNAseq enrichment.R'))
#BiocManager::install('stats')
library(stats)
### need to load or add the function first because it fails !! 
#### Prerequisites are mofa analysis 
### and rna seq enrichment 



list1=list()
list_proteins=list()
list_mirs=list()
nfactors<-MOFAobject@dimensions$K
vars_by_factor_all<-calculate_variance_explained(MOFAobject)
group=1
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]
vars_by_factor_all
mofa_enrich_rds<-paste0(outdir, '/enrichment/gse_results_mofa.Rds')

ONT='BP'

get_ranked_gene_list_mofa<-function(view, factor){
  f1<-get_weights(MOFAobject, view=view, factor=factor)
  gene_list<-f1[view][[1]] 
  gene_list
  hist(gene_list)
  
  order_ind<-order(-gene_list)
  gene_list_ord<-gene_list[order_ind,]
  names(gene_list_ord)<-rownames(gene_list)[order_ind]
  return(gene_list_ord)
}

process_mofa=TRUE
pvalueCutoff=1
nfactors=8
if (file.exists(mofa_enrich_rds)){
  list1<-loadRDS(mofa_enrich_rds)
  
}else{
  
      for (factor in 1:nfactors){
            for (view in c(  'proteomics')){
              #### Do the RNA view for whatever is high in rna 
              print(paste0(view, factor ))
              gene_list_ord<-get_ranked_gene_list_mofa(view, factor)
              gene_list_ord
              if (view=='RNA'){
               
                ### Run RNA 
                    gse_mofa <- clusterProfiler::gseGO(gene_list_ord, 
                                                       ont=ONT, 
                                                       keyType = 'ENSEMBL', 
                                                       OrgDb = 'org.Hs.eg.db', 
                                                       pvalueCutoff  = pvalueCutoff)
                  
                    list1[[factor]]<-gse_mofa
                    }
                 
          
              
                }
                  ### Run proteins 
                if (view=='proteomics'){
                  gse_protein_full <- clusterProfiler::gseGO(gene_list_ord, 
                                                             ont=ONT, 
                                                             keyType = 'SYMBOL', 
                                                             OrgDb = 'org.Hs.eg.db', 
                                                             pvalueCutoff  = pvalueCutoff)
                  list_proteins[[factor]]<-gse_protein_full
                  
                }
               #if (view=='miRNA'){
               #   mieaa_all_gsea <- rba_mieaa_enrich(test_set = gene_list_ord,
               #                                      mirna_type = "mature",
               #                                      test_type = "GSEA",
               #                                      species = 'Homo sapiens',
               #                                      sig_level=pvalueCutoff
               #   )
                #  list_mirs[[factor]]<-mieaa_all_gsea
                  
               # }
             
          
            
          

      
      }
      
      saveRDS(list(list1,list_proteins, list_mirs),mofa_enrich_rds )

      
}


list_proteins
#### Now run the prot view ? 




#gse_mofa=list1[[factor]]

## or LOAD GSE RESULTS HERE 
for (factor in 1:nfactors){
  process_mirnas=FALSE
  results_file_mofa = paste0(outdir, '/enrichment/gsego_',factor,'_')
  gse_mofa=list1[[factor]]
  write.csv(as.data.frame(gse_mofa@result), paste0(results_file_mofa, '.csv'))
  
  ### to run mofa results
  run_enrichment_plots(gse=gse_mofa, results_file = results_file_mofa)
} 


list_prote0ins[[1]]
run_ORA=FALSE
for (factor in 1:nfactors){
    results_file_mofa = paste0(outdir, '/enrichment/proteins/gsego_',factor,'_')
    dir.create(paste0(outdir, '/enrichment/proteins/'))
    gse_mofa=list_proteins[[factor]]
    gse_mofa_sig=write_filter_gse_results(gse_mofa, results_file_mofa, pvalueCutoff)
    
    
    write.csv(as.data.frame(gse_mofa@result), paste0(results_file_mofa, '.csv'))
    
    ### to run mofa results
    
    process_mirnas=FALSE
    process_mofa=TRUE
    dim(gse_mofa_sig)
    #& vars_by_factor[,'proteomics'][factor]>1
    if  (dim(gse_mofa_sig)[1]>2 ){
      print(paste(factor,'sig'))
      which(gse_mofa@result$p.adjust<0.05)
      gse_mofa
          enrich_plots<-run_enrichment_plots(gse=gse_mofa_sig,
                                             results_file=results_file_mofa, 
                                             N_DOT=15, N_EMAP = 15)    
          
    }
      
}


for (factor in 1:nfactors){
  
  results_file_mofa = paste0(outdir, '/enrichment/mirnas/gsego_',factor,'_')
  dir.create(paste0(outdir, '/enrichment/mirnas/'))
  gse_mofa=list_mirs[[factor]]
  gse_mofa_sig=write_filter_gse_results(gse_mofa, results_file_mofa, pvalueCutoff)
  
  
  write.csv(as.data.frame(gse_mofa@result), paste0(results_file_mofa, '.csv'))
  
  ### to run mofa results
  
  process_mirnas=FALSE
  process_mofa=TRUE
  dim(gse_mofa_sig)
  #& vars_by_factor[,'proteomics'][factor]>1
  if  (dim(gse_mofa_sig)[1]>2 ){
    print(paste(factor,'sig'))
    which(gse_mofa@result$p.adjust<0.05)
    gse_mofa
    enrich_plots<-run_enrichment_plots(gse=gse_mofa_sig,
                                       results_file=results_file_mofa, 
                                       N_DOT=20, N_EMAP = 15)    
    
  }
  
}


vars_by_factor[,'proteomics']
    
