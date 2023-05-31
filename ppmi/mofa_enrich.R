
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, 'ppmi/utils.R'))
## run mofa / or just load model 
source(paste0(script_dir,'ppmi/mofa_application_ppmi_all_visits.R'))
library(R.filesets)
library(stats)
### need to load or add the function first because it fails !! 
#### Prerequisites are mofa analysis 
### and rna seq enrichment 



#### Take a mofa directory and run enrichment for important factors ####
#### 1. Run for factors with cor> cor_T
#### 2. 



nfactors<-MOFAobject@dimensions$K;nfactors
vars_by_factor_all<-calculate_variance_explained(MOFAobject)
group=1
vars_by_factor<-vars_by_factor_all$r2_per_factor[[group]]
vars_by_factor_all
dir.create(paste0(outdir, '/enrichment/'))
mofa_enrich_rds<-paste0(outdir, '/enrichment/gse_results_mofa')


#cors_pearson_l<-read.csv(paste0(outdir, '/covariate_corelations_pearson.csv'))
cors_pearson_l<-correlate_factors_with_covariates(MOFAobject,
                                  covariates = c('CONCOHORT', 'INEXPAGE'), 
                                  plot = "r", 
                                  return_data = TRUE
                                  
)
cors_pval<-correlate_factors_with_covariates(MOFAobject,
                                                  covariates = c('CONCOHORT','INEXPAGE'), 
                                                  plot = "log_pval", 
                                                  return_data = TRUE
                                                  
)


cohort_cors<-cors_pearson_l[,'CONCOHORT'] # TODO: LOAD or recalc
cors_pearson_l[,'INEXPAGE'] 
cors_pval[,'CONCOHORT']
cohort_cors
cor_t=0.15
sel_factors<-which(abs(cohort_cors)>cor_t)
sel_factors
cohort_cors[sel_factors]

vars_by_factor_all$r2_per_factor$group1[sel_factors,]

ONT='BP'

get_ranked_gene_list_mofa<-function(view, factor){
  f1<-get_weights(MOFAobject, view=view, factor=factor)
  gene_list<-f1[view][[1]] 
  gene_list
  #phist<-hist(gene_list)
  #ggsave(paste0(outdir, 'pvalue_hist.jpeg'), phist)
  order_ind<-order(-gene_list)
  gene_list_ord<-gene_list[order_ind,]
  names(gene_list_ord)<-rownames(gene_list)[order_ind]
  return(gene_list_ord)
}

process_mofa=TRUE
pvalueCutoff=1
#nfactors=12
list1= vector("list", length = nfactors)
list_proteins= vector("list", length = nfactors)
list_mirs= vector("list", length = nfactors)


mofa_enrich_rds<-paste0(outdir, '/enrichment/gse_results_mofa')

#sel_factors=c(4)
#if (file.exists(mofa_enrich_rds)){
 # list_all<-loadRDS(mofa_enrich_rds)
 # list_proteins<-loadRDS(paste0(mofa_enrich_rds, 'prot'))
#  list_mirs<-loadRDS(paste0(mofa_enrich_rds, 'mirs'))
  
#}else{
### ONLY run on the server as rstudio script , there is a problem in Rstudio
#sel_factors=1:8
if (!isRStudio){
  


      for (factor in sel_factors){
           # for (view in c( 'proteomics')){
              #for (view in c( 'RNA', 'miRNA', 'proteomics')){
             #  for (view in c( 'RNA', 'miRNA', 'proteomics')){
               for (view in c( 'RNA', 'miRNA', 'proteomics')){
             #for (view in c(  'miRNA', 'proteomics')){

              #### Do the RNA view for whatever is high in rna
                    print(paste0(view,' ', factor ))
                    gene_list_ord<-get_ranked_gene_list_mofa(view, factor)
                    gene_list_ord
                    if (view=='RNA'){
                     
                      ### Run RNA 
                      if (FALSE){
                        
                         # if (file.exists(paste0(mofa_enrich_rds, 'gene'))){
                    #      if (file.exists(paste0(mofa_enrich_rds, 'gene'))){
                            ## to RERUN WITH NEW FACTORS YOU need to force it
                            list1<-loadRDS(paste0(mofa_enrich_rds, 'gene'))
                          }else{
                            gene_list_ora<-gene_list_ord[order(abs(gene_list_ord))]
                            
                                gse_mofa <- clusterProfiler::gseGO(gene_list_ord, 
                                                                   ont=ONT, 
                                                                   keyType = 'ENSEMBL', 
                                                                   OrgDb = 'org.Hs.eg.db', 
                                                                   pvalueCutoff  = pvalueCutoff)
                                
                                gse_mofa <- clusterProfiler::enrichGO(names(gene_list_ora)[1:50], 
                                                                   ont=ONT, 
                                                                   keyType = 'ENSEMBL', 
                                                                   OrgDb = 'org.Hs.eg.db', 
                                                                   pvalueCutoff  = pvalueCutoff)
                                
                                gse_mofa@result$Description[which(gse_mofa@result$p.adjust<0.05)]
                                list1[[factor]]<-gse_mofa
                                saveRDS(list1, paste0(mofa_enrich_rds, 'gene'))
                            
                          }
                    
                    
              
                }
                  ### Run proteins 
                if (view=='proteomics'){
                  
                  
                 #if (file.exists(paste0(mofa_enrich_rds, 'prot'))){
                  if (FALSE){

                    list_proteins<-loadRDS(paste0(mofa_enrich_rds, 'prot'))
                  }else{
                  
                          
                          gse_protein_full <- clusterProfiler::gseGO(gene_list_ord, 
                                                                     ont=ONT, 
                                                                     keyType = 'SYMBOL', 
                                                                     OrgDb = 'org.Hs.eg.db', 
                                                                     pvalueCutoff  = pvalueCutoff)
                          list_proteins[[factor]]<-gse_protein_full
                          saveRDS(list_proteins, paste0(mofa_enrich_rds, 'prot'))
                  
                  }
                }
              if (view=='miRNA'){
                
                
                if (FALSE){
                  
                  list_mirs<-loadRDS(paste0(mofa_enrich_rds, 'mirs'))
                
                         
                }else{
                
                        #gene_list_ord_cut<-gene_list_ord[order(abs(gene_list_ord))[1:200]]
                        #gene_list_ord_cut<-gene_list_ord_cut[order(gene_list_ord_cut, decreasing=TRUE)]
                         mieaa_all_gsea_mofa <- rba_mieaa_enrich(test_set = names(gene_list_ord),
                                                            mirna_type = "mature",
                                                            test_type = "GSEA",
                                                            species = 'Homo sapiens'                      ,
                                                            categories = c('miRPathDB_GO_Biological_process_mature'),
                                                            sig_level=pvalueCutoff
                         )
        
                         list_mirs[[factor]]<-mieaa_all_gsea_mofa
                         saveRDS(list_mirs, paste0(mofa_enrich_rds, 'mirs'))
                }
               }
             
          
            
          

      
              }
      }
}     
#      saveRDS(list(list1,list_proteins, list_mirs),mofa_enrich_rds )

      
#}
### now load already executed analysis on the server

list1<-loadRDS(paste0(mofa_enrich_rds, 'gene'))


list_mirs<-loadRDS(paste0(mofa_enrich_rds, 'mirs'))
list_proteins<-loadRDS(paste0(mofa_enrich_rds, 'prot'))

as.logical(lapply(list1, is.null))
as.logical(lapply(list_mirs, is.null))
as.logical(lapply(list_proteins, is.null))


list_mirs_enrich=list()
sel_factors_to_p<-sel_factors
sel_factors
sel_factors_to_p<-sel_factors[1:2]
sel_factors_to_p

for (factor in sel_factors){
  results_file_mofa = paste0(outdir, '/enrichment/mirnas/gsego_',factor,'_')
  suppressWarnings(dir.create(paste0(outdir, '/enrichment/mirnas/')))
  gse_mofa_mirs=list_mirs[[factor]]
  
  mieaa_all_gsea=gse_mofa_mirs
  Padj_T_paths=0.05
  mieaa_res<-mirna_enrich_res_postprocessing(mieaa_all_gsea, mir_results_file=results_file_mofa)
  mieaa_gsea_1=mieaa_res[[1]]
  enr_full=mieaa_res[[2]]
  
  
  list_mirs_enrich[[factor]]<-enr_full
}

as.logical(lapply(list_mirs_enrich, is.null))


#list1=listALL[[1]]

#### Now run the prot view ? 
run_plots=ifelse(isRStudio, FALSE, TRUE)


if (run_plots){

run_ORA=FALSE
## or LOAD GSE RESULTS HERE 


for (factor in sel_factors_to_p){
  results_file_mofa = paste0(outdir, '/enrichment/mirnas/gsego_',factor,'_')
  dir.create(paste0(outdir, '/enrichment/mirnas/'))
  gse_mofa_mirs=list_mirs[[factor]]
  
  mieaa_all_gsea=gse_mofa_mirs
  Padj_T_paths=0.05
  mieaa_res<-mirna_enrich_res_postprocessing(mieaa_all_gsea, mir_results_file=results_file_mofa)
  mieaa_gsea_1=mieaa_res[[1]]
  enr_full=mieaa_res[[2]]
  
  
  list_mirs_enrich[[factor]]<-enr_full
  
  print(any(enr_full@result$p.adjust<0.05))
  
  
  gse_mofa_sig=write_filter_gse_results(enr_full, results_file_mofa, pvalueCutoff)
  
  write.csv(as.data.frame(gse_mofa_sig@result), paste0(results_file_mofa, '.csv'))
  
  ### to run mofa results
  
  process_mirnas=TRUE
  process_mofa=TRUE
  dim(gse_mofa_sig)
  #& vars_by_factor[,'proteomics'][factor]>1
  type(gse_mofa_sig)
  dim(gse_mofa_sig)[1]>4
  if  (dim(gse_mofa_sig)[1]>4 ){
    print(dim(gse_mofa_sig)[1])
    print(paste(factor,'sig'))
    enrich_plots<-run_enrichment_plots(gse=gse_mofa_sig,
                                       results_file=results_file_mofa, 
                                       N_DOT=20, N_EMAP = 15)    
    
  }
  
}


sel_factors_to_p
for (factor in sel_factors_to_p){
  process_mirnas=FALSE
  results_file_mofa = paste0(outdir, '/enrichment/gsego_',factor,'_')
  gse_mofa_rna=list1[[factor]]
  write.csv(as.data.frame(gse_mofa_rna@result), paste0(results_file_mofa, '.csv'))
  
  ### to run mofa results
  run_enrichment_plots(gse=gse_mofa_rna, results_file = results_file_mofa)
} 


list_proteins[[1]]

process_mofa=TRUE
for (factor in sel_factors_to_p){
    results_file_mofa = paste0(outdir, '/enrichment/proteins/gsego_',factor,'_')
    dir.create(paste0(outdir, '/enrichment/proteins/'))
    gse_mofa=list_proteins[[factor]]
    gse_mofa_sig=write_filter_gse_results(gse_mofa, results_file_mofa, 0.2)
    
    
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




}

