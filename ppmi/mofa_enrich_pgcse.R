
library('MOFAdata')
utils::data(reactomeGS)

head((reactomeGS))

## TODO: if enrichment is already run then just load results
## load res.positive to be used in the next script
res.positive$feature.sets
res=res.positive
res.positive
write_enrich<-function(res, sign_mode){
  #' 
  #'' @res res.negative result from mofa enrichment 
  #'
  #'
  #'
      results_enrich<-res$pval.adj
      all_fs_merged2<-reshape::melt(results_enrich)
      all_fs_merged2_pval<-reshape::melt(res$pval )
      all_fs_merged2_pval2<-merge(all_fs_merged2,all_fs_merged2_pval , by=c('X1', 'X2'))
      #all_fs_merged2<-all_fs_merged2[all_fs_merged2$value<T,]
      all_fs_merged2_pval2<-all_fs_merged2_pval2[with(all_fs_merged2_pval2, order(X2, value.x)),]# order 
      colnames(all_fs_merged2_pval2)<-c('Description', 'Factor', 'p.adjust', 'pvalue')
      
      neg_file<-paste0(outdir,'/enrichment/',gsub('\\:', '_', subcategory), 
                       mode, '_enrichment', sign_mode)
      write.csv(format(all_fs_merged2_pval2, digits=3),paste0(neg_file,  '.csv' ))
      all_fs_merged2_pval2_sig=all_fs_merged2_pval2[ all_fs_merged2_pval2$p.adjust<T,]
      write.csv(format(all_fs_merged2_pval2_sig, digits=3),paste0(neg_file, '_', T,  '.csv' ))
      return(all_fs_merged2_pval2)
}


subcategory<- 'CP:KEGG'
subcategory<- 'CP:KEGG'
subcategory<- 'GO:MF'
subcategory<- 'GO:BP'
dir.create(paste0(outdir, '/enrichment/'))
#for (subcategory in c('GO:BP' ,'CP:KEGG')){

mode='proteomics'
mode='RNA'

for (subcategory in c('GO:BP', 'GO:MF' )){
  if (mode=='proteomics'){
    gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), 'proteins.csv')
    
  }else{
    gs_file<-paste0(output_files, 'gs', gsub('\\:', '_', subcategory), '.csv')
    
  }
  
  gs<-as.matrix(read.csv(gs_file, header=1, row.names=1))
  colnames(gs)
  
  
  features_names(MOFAobject)$RNA
  features_names(MOFAobject)$RNA<-sapply(features_names(MOFAobject)$RNA, 
                                         function(x) {stringr::str_remove(x, '\\..*')}
  )
  
  
  sign_mode='negative'
  enrich_res_file_neg<-paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'negative' )
  enrich_res_file_pos<-paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'positive' )
  
  if (file.exists(enrich_res_file_neg)){
        res.negative=loadRDS(enrich_res_file_neg)
        res.positive=loadRDS(enrich_res_file_pos)
    
  }else{
          # GSEA on negative weights, with default options
          res.negative <- run_enrichment(MOFAobject, 
                                         feature.sets = gs, 
                                         view = mode,
                                         sign = "negative"
          )
          
          res.positive <- run_enrichment(MOFAobject, 
                                         feature.sets = gs, 
                                         view = mode,
                                         sign = "positive"
          )
          
          
          res_negative_df<-write_enrich(res.negative, sign_mode='negative')
          saveRDS(res.negative, paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'negative' ))
          
          res_positive_df<-write_enrich(res.positive, sign_mode='positive')
          saveRDS(res.positive, paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment_', 'positive' ))
          res_merged<-merge(res_negative_df, res_positive_df, suffixes=c('_n', '_p'),
            by=c('Description','Factor' ))
          res_merged$pvalue_min<-c(rowMins(as.matrix(res_merged[,c('pvalue_n','pvalue_p')])))
          
          write.csv(res_merged, paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment.csv' ), row.names=FALSE)
  }
  
  
  
  
  

    ## TODO: create a function to do for both positive and negative 
    #
    T=0.05
  
 
  
  
  
  
  
  
  ##### which factor is related to parkinsons disease in KEGG
  ### PROBLEM: this is based on RNA only!!! 
  
  
}
sign_mode='negative'
subcategory<- 'GO:BP'
T=0.05

res.negative_l$pval

# Make enrichment plots for all factors 
# threshold on p value to zoom in 
jpeg(paste0(outdir,'/enrichment/Enrichment_heatmap_positive','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.positive, 
                        alpha=0.5, 
                        cap=0.0005,
                        colnames=TRUE)
dev.off()

plot_enrichment_heatmap(res.positive$sigPathways, 
                        alpha=0.5, cap=0.0005)

#ggsave(paste0(outdir,'Enrichment_heatmap_positive','.jpeg'), width = 9, height=4, dpi=120)


jpeg(paste0(outdir,'/enrichment/Enrichment_heatmap_negative','.jpeg'), res=150, height=800, width=800)

plot_enrichment_heatmap(res.negative, 
                        alpha=0.5, cap=0.00000000005 
)


dev.off()
#ggsave(paste0(outdir,'Enrichment_heatmap_negative','.png'), width = 9, height=4, dpi=120)


F3<-res.positive$pval.adj[,'Factor3']
SIG<-F3[F3<0.05]
SIG[order(SIG)][1:20]

F3<-res.negative$pval.adj[,'Factor6']
SIG<-F3[F3<0.05]
SIG[order(SIG)][1:10]
