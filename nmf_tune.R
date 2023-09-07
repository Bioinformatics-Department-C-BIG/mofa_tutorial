


unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

unregister_dopar()

run_nmf_get_cors<-function(){
  
  output_cors=list()
  #for (g_params in c('0.05_100_', '0.1_100_', '0.2_100_', '0.3_100_')){
  #  for (g_params in c( '0.3_100_')){
  for (g_params in c( '0.4_100_')){
    
    for (NFACTORS in seq(5,11)){
      
      for  (k_centers in c(3,4,5,6,7,8,9)){
         print(k_centers)
      

            
            print(paste(g_params, NFACTORS))
            #### Load the dataset ####
            ## get list of three mats 
            param_str_g<-paste0('rnas_', VISIT_S, '_', g_params, 'coh_', sel_coh_s, '_'  , sel_subcoh_s)
            
            data_full=prepare_multi_data(p_params, param_str_g, param_str_m, mofa_params)
            # create multi experiment 
            mofa_multi<-create_multi_experiment(data_full, combined_bl)
            
            mofa_multi$COHORT
            
            
            x1_se<-mofa_multi[, , mod]
            
            x1=assays(x1_se)[[mod]]
            
            print(NFACTORS)
            nrun=3
            
            
            ## SAVE AND LOAD 
           # out_nmf<-paste0('nmf/plots/','multirun_', param_str_g, NFACTORS, '_', nrun,'_', 
           #                 mod)
            out_nmf_params<- paste0( 'g_', g_params,'_coh_', sel_coh_s,'_', VISIT_S)
            outdir_nmf = paste0(outdir_orig, '/nmf/',out_nmf_params, '_', NFACTORS );
            out_nmf=paste0(outdir_nmf, 'model')
            seed=135
            
            if (file.exists(out_nmf)){
              res=loadRDS(out_nmf)
              
            }else{
              
              res.multirun<-NMF::nmf(x1,NFACTORS,nrun=nrun, seed=seed )
              res=res.multirun
              saveRDS(res.multirun,out_nmf)
            }
            
            
            ### return fitted model 
            fit(res)
            h <-as.data.frame(coef(res)) # factor coeficients for each sample 
            dim(h)
            h[1,]
            x1_se$PATNO
            
            #### Save output correlations 

            cor_res<-round(apply(h, 1, cor, x=x1_se$COHORT), 2)
            cor_res
           # output_cors[[NFACTORS]]<-cor_res
            
            
            
            
            ### Also cluster and measure cluster corelation ###
            
            
            h <-as.data.frame(coef(res)) # factor coeficients for each sample 
    
            #### Correlations for each factor
            # 1. Tune to maximize cohort correlations ####
            covariates <- as.data.frame(lapply(colData(x1_se), as.numeric))
            rownames(covariates)<-colData(x1_se)$PATNO_EVENT_ID
            h_t<-as.data.frame(t(h))
            
            
            
            #### 
            cor <- psych::corr.test(covariates,h_t, method = "pearson", adjust = "BH")
            print(cor$p.adj['COHORT',]<0.05)
            
            
            
            
            
            
            
            
            sel_factors<-which(cor$p['COHORT',]<0.05)
            
          
            cor_clusters=NULL; mut_inf=NULL
           
            print(sel_factors)
            print(paste0('kcenters: ', k_centers))
            if (length(sel_factors)>0){
              
              clusters_single <- kmeans(t(h)[,sel_factors], centers = k_centers)
              
              covariates$cluster_s<-clusters_single$cluster[match(rownames(covariates),names(clusters_single$cluster))]

              
              
              cor_clusters<-chisq.test(clusters_single$cluster,covariates$COHORT )
              print(cor_clusters)
              cor_clusters$p.value
              
              mut_inf<-round( MutInf(clusters_single$cluster,covariates$COHORT),digits = 3)
              print(paste('MI:', mut_inf, digits = 4))
              
            }
              
            #df_stats=  c( TOP_PN, TOP_GN, MIN_COUNT_G, TOP_MN, MIN_COUNT_M, mofa_params, sel_coh_s,VISIT_S,  scale_views[1],  use_signif,
            #              run_mofa_complete, N_FACTORS,cors_t , max_cor )
            
            
            df_stats = c(g_params, NFACTORS, k_centers, cor_clusters$p.value, cor_clusters$statistic, mut_inf, length(sel_factors), nrun )
            write.table(t(df_stats), paste0(outdir_orig,'nmf_all_stats_mi.csv'), append=TRUE,sep=',', col.names = FALSE)
            
          }

  
  
      }
  }
}



run_nmf_get_cors()

