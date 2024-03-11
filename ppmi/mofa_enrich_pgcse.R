

#library('MOFAdata')
library('MOFA2')



## TODO: if enrichment is already run then just load results
## load res.positive to be used in the next script



write_enrich<-function(res, sign_mode){
  #' 
  #'' @res res.negative result from mofa enrichment 
  #'
  #'
  #'res$pval.adj
      res$pval.adj=as.data.frame(res$pval.adj)   ;
      results_enrich<-as.data.frame(sapply(res$pval.adj, as.numeric));
      results_enrich$Description=rownames(res$pval.adj); 
      all_fs_merged2<-reshape::melt(results_enrich, id.vars = 'Description');
      
      res$pval=as.data.frame(res$pval) 
      res_enrich_pval=as.data.frame(sapply(res$pval, as.numeric))
      res_enrich_pval$Description=rownames(res$pval);
      
      all_fs_merged2_pval<-reshape::melt(res_enrich_pval, id.vars = 'Description');
      all_fs_merged2_pval2=cbind(all_fs_merged2, all_fs_merged2_pval$value)
      #rownames(all_fs_merged2_pval2)=rownames(all_fs_merged2_pval);
      
      #all_fs_merged2_pval2<-merge(all_fs_merged2,all_fs_merged2_pval , by=c('variable', 'value'))

      rownames(all_fs_merged2_pval2)
      all_fs_merged2_pval2<-all_fs_merged2_pval2[with(all_fs_merged2_pval2, order(variable, value)),]# order 
     colnames(all_fs_merged2_pval2)<-c('Description', 'Factor', 'p.adjust', 'pvalue')
      #colnames(all_fs_merged2_pval2)<-c('Factor', 'p.adjust', 'pvalue')
      
      neg_file<-paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory), 
                       mode, '_enrichment', sign_mode)
      
      write.csv(format(all_fs_merged2_pval2, digits=3),paste0(neg_file,  '.csv' ), row.names = TRUE)
      all_fs_merged2_pval2_sig=all_fs_merged2_pval2[ all_fs_merged2_pval2$p.adjust<T,]
      write.csv(format(all_fs_merged2_pval2_sig, digits=3),paste0(neg_file, '_', T,  '.csv' ))
      
}


subcategory<- 'CP:KEGG'
subcategory<- 'CP:KEGG'
subcategory<- 'GO:MF'
subcategory<- 'GO:BP'
dir.create(paste0(outdir, '/enrichment/'))
#for (subcategory in c('GO:BP' ,'CP:KEGG')){

mode='proteomics'
mode='RNA'

mode='proteomics_csf'
features_names(MOFAobject)$proteomics_csf


features_names(MOFAobject)$proteomics_csf
modify_feature_names<-function(MOFAobject){

    #original=features_names(MOFAobject)$proteomics_csf
    features_names(MOFAobject)$proteomics_csf<-gsub('_proteomics_csf', '', features_names(MOFAobject)$proteomics_csf)
#      features_names(MOFAobject)$proteomics_csf<-gsub('_c', '', features_names(MOFAobject)$proteomics_csf)

    #features_names(MOFAobject)$proteomics_csf<-gsub('_c', '', rownames(MOFAobject@expectations$W$proteomics_csf))
    features_names(MOFAobject)$proteomics_plasma<-gsub('_proteomics_plasma', '_', features_names(MOFAobject)$proteomics_plasma)
   return(MOFAobject)

}


# rename to match databases 
MOFAobject_enr<-MOFAobject
MOFAobject_enr<-modify_feature_names(MOFAobject)
head(features_names(MOFAobject_enr)$proteomics_csf)
head(features_names(MOFAobject_enr)$proteomics_plasma)

head(features_names(MOFAobject_enr)$RNA)
features_names(MOFAobject_enr)$RNA<-sapply(features_names(MOFAobject_enr)$RNA, 
                                         function(x) {stringr::str_remove(x, '\\..*')}
  )



for (subcategory in c('GO:BP', 'GO:MF' )){
  if (mode=='proteomics_csf'| mode=='proteomics_plasma'){
    gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), 'proteins.csv')
    
  }else{
    gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
    
  }
  
  gs<-as.matrix(read.csv(gs_file, header=1, row.names=1))
  colnames(gs)
  
  

  
  
  sign_mode='negative'
  enrich_res_file_neg<-paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'negative' )
  enrich_res_file_pos<-paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'positive' )
  
  if (file.exists(enrich_res_file_neg)){
        res.negative=loadRDS(enrich_res_file_neg)
        res.positive=loadRDS(enrich_res_file_pos)
    
  }else{
          # GSEA on negative weights, with default options
          res.negative <- run_enrichment(MOFAobject_enr, 
                                         feature.sets = gs, 
                                         view = mode,
                                         sign = "negative"
          )
          
          res.positive <- run_enrichment(MOFAobject_enr, 
                                         feature.sets = gs, 
                                         view = mode,
                                         sign = "positive"
          )
          
          sign_mode='negative'
          res_negative_df<-write_enrich(res.negative, sign_mode=sign_mode)
          saveRDS(res.negative, paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'negative' ))
          

          sign_mode='positive'

          res_positive_df<-write_enrich(res.positive, sign_mode=sign_mode)
          saveRDS(res.positive, paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'positive' ))
          res_merged<-merge(res_negative_df, res_positive_df, suffixes=c('_n', '_p'),
            by=c('Description','Factor' ))
          res_merged$pvalue_min<-c(rowMins(as.matrix(res_merged[,c('pvalue_n','pvalue_p')])))
          
          write.csv(res_merged, paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment.csv' ), row.names=FALSE)
  }
  
 
  
  
  
}



























