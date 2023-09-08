


### Complex heatmap 

library('ComplexHeatmap')


plot_heatmap<-function(vsd_filt, sigGenes,  df,remove_cn=FALSE, show_rownames=TRUE, 
                       cluster_cols=FALSE, order_by_hm='COHORT'){
  #'
  #' @param vsd_filt: annotation dataframe nsamples X ncoldata 
  #' @param hm: heatmap data nfeats X nsamples 
  #' @param sigGenes: select genes to plot description
  #' 
  #' 
  #' 
  #' 
  #ARRANGE
  
  
  
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
  
  
  ### Add the annotations 
  meta_single<-colData(vsd_filt_genes)
  
  
  
  
  ## HEATMAP OPTIONS 
  cluster_cols=cluster_cols;    
  #colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
  fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
                filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,remove_cn ,'.jpeg')
  
  hm<-assay(vsd_filt_genes)
  
  df_ord<-df[order(df[,order_by_hm]),]
  df_ord$COHORT
  dim(hm)
  hm_ord<-hm[,order(df[,order_by_hm])]
  
  ### SCALE!! 
  hm_scaled <- as.matrix(hm_ord) %>% t() %>% scale() %>% t()
  dim(hm_ord)
  cluster_cols=cluster_cols
  hm_scaled
  #jpeg(fname, width=2000, height=1500, res=200)
  graphics.off()
  library(ggplot2)
  if(process_mirnas){
    lab=rownames(rowData(vsd_filt_genes)) }else{
      lab=as.character(rowData(vsd_filt_genes)$SYMBOL)}
  cluster_cols
  hm_scaled
  
  hm_scaled_filt=hm_scaled
  df_ord_filt=df_ord
  
  if (remove_cn){
    d_ind<-df_ord$COHORT==1
    hm_scaled_filt<-hm_scaled[,d_ind]
    
    df_ord_filt<-df_ord[d_ind, ]
  }
  
  #jpeg(fname, width=10*100, height=7*100, res=300)
  my_pheatmap<-pheatmap::pheatmap(hm_scaled_filt, 
                        labels_row=lab,
                        cluster_rows=TRUE, 
                        show_rownames=show_rownames,
                        #scale='row', 
                        cluster_cols=cluster_cols,
                        
                        annotation_col=df_ord_filt, 
                        
                        clustering_method = 'ward.D2'
  )
  
  
  show(my_pheatmap)
  dev.off()
  my_pheatmap
  # ggsave(fname, width=10, height=7)
  
  ggsave(fname,plot=my_pheatmap, width=7, height=7, dpi=300)
  return(my_pheatmap)
}




remove_cn=FALSE
order_by_hm='COHORT'

plot_heatmap<-function(vsd_filt, sigGenes,  df,remove_cn=FALSE, show_rownames=TRUE, 
                       cluster_cols=FALSE, order_by_hm='COHORT'){
  
  
  
  #'
  #' @param vsd_filt: annotation dataframe nsamples X ncoldata 
  #' @param hm: heatmap data nfeats X nsamples 
  #' @param sigGenes: select genes to plot description
  #' 
  #' 
  #' 
  #' 
  #ARRANGE
  
  
  
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
  
  
  ### Add the annotations 
  meta_single<-colData(vsd_filt_genes)
  
  
  
  
  ## HEATMAP OPTIONS 
  cluster_cols=cluster_cols;   
  
  #colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
  fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
                filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,remove_cn ,'.jpeg')
  
  hm<-assay(vsd_filt_genes)
  
  df_ord<-df[order(df[,order_by_hm]),]
  df_ord$COHORT
  dim(hm)
  hm_ord<-hm[,order(df[,order_by_hm])]
  
  ### SCALE!! 
  hm_scaled <- as.matrix(hm_ord) %>% t() %>% scale() %>% t()
  dim(hm_ord)
  cluster_cols=cluster_cols
  hm_scaled
  #jpeg(fname, width=2000, height=1500, res=200)
  library(ggplot2)
  if(process_mirnas){
    lab=rownames(rowData(vsd_filt_genes)) }else{
      lab=as.character(rowData(vsd_filt_genes)$SYMBOL)}
  cluster_cols
  hm_scaled
  
  hm_scaled_filt=hm_scaled
  df_ord_filt=df_ord
  
  if (remove_cn){
    d_ind<-df_ord$COHORT==1
    hm_scaled_filt<-hm_scaled[,d_ind]
    
    df_ord_filt<-df_ord[d_ind, ]
  }
  
  #jpeg(fname, width=10*100, height=7*100, res=300)
  my_pheatmap<-pheatmap(hm_scaled_filt, 
                        labels_row=lab,
                        cluster_rows=TRUE, 
                        show_rownames=show_rownames,
                        #scale='row', 
                        cluster_cols=cluster_cols,
                        
                        annotation_col=df_ord_filt, 
                        
                        clustering_method = 'ward.D2'
  )
  
  
  
  
  ComplexHeatmap::pheatmap(hm_scaled_filt, 
                           labels_row=lab,
                           cluster_rows=TRUE, 
                           show_rownames=show_rownames,
                           #scale='row', 
                           cluster_cols=cluster_cols,
                           
                           annotation_col=df_ord_filt, 
                           
                           clustering_method = 'ward.D2')
  
  
  dev.off()
  # ggsave(fname, width=10, height=7)
  
  ggsave(fname,plot=my_pheatmap, width=7, height=7, dpi=300)
  return(my_pheatmap)
}



