
source(paste0('ppmi/setup_os.R'))
source(paste0(script_dir, 'ppmi/utils.R'))
## run mofa / or just load model 
#cell_corr=TRUE

 #source(paste0(script_dir,'ppmi/mofa_application_ppmi_all_visits.R'))
#source(paste0(script_dir,'ppmi/mofa_analysis_time_diff.R'))

library(R.filesets)
library(stats)
### need to load or add the function first because it fails !! 
#### Prerequisites are mofa analysis 
### and rna seq enrichment 



#### Take a mofa directory and run enrichment for important factors ####
#### 1. Run for factors with cor> cor_T
#### 2. 


### Takes outdir and MOFA object as input
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

cors<-cors_pval[,'CONCOHORT'] # TODO: LOAD or recalc
cohort_cors<-cors_pearson_l[,'CONCOHORT'] # TODO: LOAD or recalc
cors_pearson_l[,'INEXPAGE'] 
cor_t=0.15
sel_factors<-which(abs(cohort_cors)>cor_t)
sel_factors
cohort_cors[sel_factors]
## select by correlation? 
sel_factors<-which(cors>-log10(0.05))
sel_factors

vars_by_factor_all$r2_per_factor$group1[sel_factors,]

ONT='BP'
# todo: add ONT in results file settings 

get_ranked_gene_list_mofa<-function(view, factor){
        f1<-get_weights(MOFAobject, view=view, factor=factor)
        gene_list<-f1[view][[1]] 
        gene_list
        #phist<-hist(gene_list)
        #ggsave(paste0(outdir, 'pvalue_hist.jpeg'), phist)
        order_ind<-order(-gene_list)
        gene_list_ord<-gene_list[order_ind,]
        names(gene_list_ord)<-rownames(gene_list)[order_ind]
        names(gene_list_ord)<-gsub('\\..*', '', names(gene_list_ord)) 
        
        return(gene_list_ord)
}

process_mofa=TRUE
pvalueCutoff=1
#nfactors=12
list1= vector("list", length = nfactors)
list_proteins= vector("list", length = nfactors)
list_proteins_enrich= vector("list", length = nfactors)

list_mirs= vector("list", length = nfactors)
list_mirs_targets= vector("list", length = nfactors)

list1_genes= vector("list", length = nfactors)
list_proteins_weights=vector("list", length = nfactors)

suppressWarnings(dir.create(paste0(outdir, '/enrichment/')))

mofa_enrich_rds<-paste0(outdir, '/enrichment/gse_results_mofa')

#sel_factors=c(4)
#if (file.exists(mofa_enrich_rds)){
 # list_all<-loadRDS(mofa_enrich_rds)
 # list_proteins<-loadRDS(paste0(mofa_enrich_rds, 'prot'))
#  list_mirs<-loadRDS(paste0(mofa_enrich_rds, 'mirs'))
  
#}else{
### ONLY run on the server as rstudio script , there is a problem in Rstudio
sel_factors_to_enrich<-sel_factors
sel_factors_to_enrich=1:15
sel_factors_to_enrich<-sel_factors
sel_factors_to_enrich=1:N_FACTORS
#BiocManager::install("fgsea")


just_load=FALSE
just_load=FALSE

if (TRUE){
  
      for (factor in sel_factors_to_enrich){
          for (view in c('RNA', 'proteomics_t_csf', 'proteomics_t_plasma')){
          #view='RNA'; factor=3
                    print(paste0(view,' ', factor ))
                    #factor=4;view='proteomics'
                    gene_list_ord<-get_ranked_gene_list_mofa(view, factor)
                    if (view=='RNA'){
                     
                      ### Run RNA 
                      if (just_load){
                  
                            list1<-loadRDS(paste0(mofa_enrich_rds, 'gene'))
                            list1_genes[[factor]]<-gene_list_ord
                            
                          }else{
                          gene_list_ora<-gene_list_ord[order(abs(gene_list_ord))]
                          results_file_mofa = paste0(outdir, '/enrichment/rnas_/gsego_', ONT, '_',factor,'_' )
                          
                          gse_mofa<-run_enrich_gene_list(gene_list_ord, results_file=results_file_mofa,N_NET=20  )
                               
                              
                                list1[[factor]]<-gse_mofa
                                list1_genes[[factor]]<-gene_list_ord
                                saveRDS(list1, paste0(mofa_enrich_rds, 'gene'))
                            
                          }
                    
                    
              
                }
                  ### Run proteins 
                if (view=='proteomics_t_csf' & just_load){
                  
                 #if (file.exists(paste0(mofa_enrich_rds, 'prot'))){
                    list_proteins<-loadRDS(paste0(mofa_enrich_rds, 'prot_t_csf'))
                    tissue = 'proteomics_t_csf'
                    
                }else if  (view=='proteomics_t_plasma' & just_load){
                    list_proteins_plasma<-loadRDS(paste0(mofa_enrich_rds, 'prot_t_plasma'))
                    tissue = 'proteins_t_plasma'
                    
                  }
                    
                    
                    ### set filenames 
       
                    if (view=='proteomics_t_csf' ){
                      fname<-paste0(mofa_enrich_rds, 'prot_csf')
                      tissue = 'proteins_csf'
                      
                    }else if(view=='proteomics_t_plasma'){
                      fname<-paste0(mofa_enrich_rds, 'prot_t_plasma')
                      tissue = 'proteins_t_plasma'
                      
                                          }
                    
                    
                    
                    if (view=='proteomics_csf'|| view=='proteomics_plasma' ){
                    
                    
                        run_ORA=TRUE
                         # TODO: DECIDE on the number
                          if (run_ORA){
                            names(gene_list_ord) = gsub('_proteomics.*', '',names(gene_list_ord))
                            results_file_mofa = paste0(outdir, '/enrichment/', tissue,'/ora/gsego_',ONT, factor,'_')
                            dir.create(paste0(outdir, '/enrichment/', tissue,'/ora/'), recursive = TRUE)
                            
                            
                            gse_protein_full_enrich =   run_ora_gene_list(gene_list_ord,results_file=results_file_mofa)
                          
                            list_proteins_enrich[[factor]]<-gse_protein_full_enrich
                            saveRDS(list_proteins_enrich, paste0(mofa_enrich_rds, 'prot_enrich_go'))
                          
                                                    
                           
                          
                          }else{
                            
                          list_proteins_weights[[factor]]<-gene_list_ord
                            
                          names(gene_list_ord) = gsub('_proteomics.*', '',names(gene_list_ord))
                           gse_protein_full <- clusterProfiler::gseGO(gene_list_ord, 
                                                                     ont=ONT, 
                                                                     keyType = 'SYMBOL', 
                                                                  OrgDb = 'org.Hs.eg.db', 
                                                                    pvalueCutoff  = pvalueCutoff)
                          
                        
                           list_proteins[[factor]]<-gse_protein_full
                          saveRDS(list_proteins,fname)
                          }
                  }
                  }
                   
      }
  

  gene_list_ord
  just_load=FALSE
  sel_factors_to_enrich
sel_factors_to_enrich = c(1:N_FACTORS)
  ######## MIRNAS ON THEIR OWN
for (factor in sel_factors_to_enrich){
          for (view in c('miRNA')){
            #view='RNA'; factor=3
            print(paste0(view,' ', factor ))
            #factor=4;view='proteomics'
            gene_list_ord<-get_ranked_gene_list_mofa(view, factor)   
              if (view=='miRNA'){
                if (just_load){
                  list_mirs<-loadRDS(paste0(mofa_enrich_rds, 'mirs'))
                }else{
              
                    gene_list_ord_cut<-gene_list_ord[order(abs(gene_list_ord))[1:200]]
                    gene_list_ord_cut<-gene_list_ord_cut[order(gene_list_ord_cut, decreasing=TRUE)]
                    mieaa_all_gsea_mofa<-run_enrich_mirnas(gene_list_ord, pvalueCutoff = pvalueCutoff, test_type='GSEA')
                    list_mirs[[factor]]<-mieaa_all_gsea_mofa
                    #list_mirs_target[[factor]]<-mieaa_all_gsea_mofa
            
                    saveRDS(list_mirs, paste0(mofa_enrich_rds, 'mirs'))
                    #saveRDS(list_mirs_target, paste0(mofa_enrich_rds, 'mirs_targets'))
                         
                }
               }
          }
      }
}
#  
list_mirs
list_mirs<-loadRDS(paste0(mofa_enrich_rds, 'mirs'))

length(list_mirs)
is.null(list_mirs[[3]])
as.logical(lapply(list1, is.null))
as.logical(lapply(list_mirs, is.null))

## holds the pvalues as an enrich result
list_mirs_enrich=list()
sel_factors_to_p<-sel_factors_to_enrich
#sel_factors_to_p<-sel_factors

dir.create(paste0(outdir, '/enrichment/'))



as.logical(lapply(list_mirs_enrich, is.null))
## this will be used further-make sure it is produced correctly
list_all=list(list1,list_proteins, list_mirs_enrich)

#list1=listALL[[1]]

#### Now run the prot view ? 
# TODO: PASS 
run_plots=TRUE
run_plots=TRUE
sel_factors_to_p
                
process_mirnas=TRUE
process_mofa=TRUE

if (run_plots){

          run_ORA=TRUE
          ## or LOAD GSE RESULTS HERE 
          
          list1_genes
          
          for (factor in sel_factors_to_p){
            results_file_mofa = paste0(outdir, '/enrichment/mirnas_/gsego_',ONT,'_' ,factor,'_')
            suppressWarnings(dir.create(paste0(outdir, '/enrichment/mirnas_/')))
            gse_mofa_mirs=list_mirs[[factor]]
            
            mieaa_all_gsea=gse_mofa_mirs
            Padj_T_paths=0.05

            if (!is.null(mieaa_all_gsea)){
                print(paste('plotting, ',factor))
                mieaa_res<-mirna_enrich_res_postprocessing(mieaa_all_gsea, mir_results_file=results_file_mofa)
                #mieaa_gsea_1=mieaa_res[[1]] # gsea result # not used anywhere...
                enr_full=mieaa_res[[2]] 
                              
                list_mirs_enrich[[factor]]<-enr_full
                gse_mofa_sig=write_filter_gse_results(enr_full, results_file_mofa, pvalueCutoff)
                
                #write.csv(as.data.frame(gse_mofa_sig@result), paste0(results_file_mofa, '.csv'))
  
                if  (dim(gse_mofa_sig)[1]>4 ){
                  print(dim(gse_mofa_sig)[1])
                  print(paste(factor,'sig'))
                  enrich_plots<-run_enrichment_plots(gse=gse_mofa_sig,
                                                    results_file=results_file_mofa, 
                                                    N_DOT=20, N_EMAP = 50)    
                  
                }
            }
            
          }
          
          
      
   



}




all_target_scan<-mget(names(mirs), revmap(targetscan.Hs.egTARGETS))
tars<-all_target_scan[[1]]

