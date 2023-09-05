




run_nmf_get_cors<-function(){
  
  output_cors=list()
  
  for (NFACTORS in seq(3,7)){
    print(NFACTORS)
    nrun=5
    
    
    ## SAVE AND LOAD 
    out_nmf<-paste0(outdir,'/../','multirun_', NFACTORS, '_', nrun,'_', 
                    mod)
    
    
    
    
    if (file.exists(out_nmf)){
      res=loadRDS(out_nmf)
      
      
    }else{
      
      res.multirun<-NMF::nmf(x1,NFACTORS,nrun=nrun )
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
    
    print(round(apply(h, 1, cor, x=x1_se$COHORT), 2))
    cor_res<-round(apply(h, 1, cor, x=x1_se$COHORT), 2)
    output_cors[[NFACTORS]]<-cor_res
    
  }

  
  
  
}

run_nmf_get_cors()
