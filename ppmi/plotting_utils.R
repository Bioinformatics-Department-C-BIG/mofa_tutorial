


### Complex heatmap 

library('ComplexHeatmap')




plot_heatmap<-function(vsd_filt, sigGenes,  df,remove_cn=FALSE, show_rownames=TRUE, 
                       cluster_cols=FALSE, order_by_hm='COHORT', sel_samples){
  
  
  
  #'
  #' @param vsd_filt: annotation dataframe nsamples X ncoldata 
  #' @param hm: heatmap data nfeats X nsamples 
  #' @param sigGenes: select genes to plot description
  #' 
  #' 
  #' 
  #' 
  #ARRANGE
  sigGenes = make.names(sigGenes)
  rownames(vsd_filt) = make.names(rownames(vsd_filt))
  
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
  
  if (remove_cn){
    d_ind<-df_ord$COHORT==1
    hm_scaled_filt<-hm_scaled[,d_ind]
    df_ord_filt<-df_ord[d_ind, ]
    
    
    hm_scaled=hm_scaled_filt
    df_ord=df_ord_filt
  }
  
  
  if (!is.null(sel_samples)){
    df_pats<-gsub('\\_.*', '',rownames(df_ord) )
    d_ind<-df_pats%in% sel_samples
    hm_scaled_filt<-hm_scaled[,d_ind]
    df_ord_filt<-df_ord[d_ind, ]
    
    hm_scaled=hm_scaled_filt
    df_ord=df_ord_filt
  }
  dim(hm_scaled_filt)
  dim(df_ord)
  show_rownames=TRUE
  #jpeg(fname, width=10*100, height=7*100, res=300)
  my_pheatmap<-pheatmap(hm_scaled, 
                        labels_row=lab,
                        cluster_rows=TRUE, 
                        show_rownames=show_rownames,
                        #scale='row', 
                        cluster_cols=cluster_cols,
                        
                        annotation_col=df_ord, 
                        
                        clustering_method = 'ward.D2'
  )
  
  
  my_pheatmap
  ggsave(fname, width=7, height=7, dpi=300)
  
  ComplexHeatmap::pheatmap(hm_scaled, 
                           labels_row=lab,
                           cluster_rows=TRUE, 
                           show_rownames=show_rownames,
                           #scale='row', 
                           cluster_cols=cluster_cols,
                           
                           annotation_col=df_ord, 
                           
                           clustering_method = 'ward.D2')
  
  
  # ggsave(fname, width=10, height=7)
  
  dev.off()
  
  return(my_pheatmap)
}







plot_heatmap_time<-function(vsd_filt, sigGenes,  df,remove_cn=FALSE, show_rownames=TRUE, 
                       cluster_cols=FALSE, order_by_hm='COHORT', sel_samples){
  
  
  
  #'
  #' @param vsd_filt: annotation dataframe nsamples X ncoldata 
  #' @param hm: heatmap data nfeats X nsamples 
  #' @param sigGenes: select genes to plot description
  #' 
  #' 
  #' 
  #' 
  #ARRANGE
  sigGenes = make.names(sigGenes)
  rownames(vsd_filt) = make.names(rownames(vsd_filt))
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
  
  
  ### Add the annotations 
  meta_single<-colData(vsd_filt_genes)
  
  
  
  
  ## HEATMAP OPTIONS 
  cluster_cols=cluster_cols;   
  
  #colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
  fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
                filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,remove_cn ,'.jpeg')
  
  hm<-assay(vsd_filt_genes); dim(hm);dim(df)

  
  new_ord<-order(df[,'EVENT_ID'], df[,'PATNO'])
  df_ord<-df[new_ord ,]
      
  df_ord$COHORT
  dim(hm); dim(df_ord)
  hm_ord<-hm[,new_ord]
  
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
  
  if (remove_cn){
    d_ind<-df_ord$COHORT==1
    hm_scaled_filt<-hm_scaled[,d_ind]
    df_ord_filt<-df_ord[d_ind, ]
    
    
    hm_scaled=hm_scaled_filt
    df_ord=df_ord_filt
  }
  
  
  if (!is.null(sel_samples)){
    df_pats<-gsub('\\_.*', '',rownames(df_ord) )
    d_ind<-df_pats%in% sel_samples
    hm_scaled_filt<-hm_scaled[,d_ind]
    df_ord_filt<-df_ord[d_ind, ]
    
    hm_scaled=hm_scaled_filt
    df_ord=df_ord_filt
  }
  dim(hm_scaled_filt)
  dim(df_ord)
  show_rownames=TRUE
  #jpeg(fname, width=10*100, height=7*100, res=300)
  
  ## Put a blank row in each one 
  split_time<-as.numeric(cumsum(table(df_ord$EVENT_ID)))
  
  df_ord<-df_ord[,!colnames(df_ord) %in% c('PATNO_EVENT_ID', 'PATNO')]
  
  
  my_pheatmap<-pheatmap(hm_scaled, 
                        labels_row=lab,
                        cluster_rows=TRUE, 
                        show_rownames=show_rownames,
                        gaps_col =  split_time,
                        #scale='row', 
                        cluster_cols=FALSE,
                        
                        annotation_col=df_ord, 
                        
                        clustering_method = 'ward.D2'
  )
  
  
  my_pheatmap
  ggsave(fname, width=7, height=7, dpi=300)
  
  
  
  ComplexHeatmap::pheatmap(hm_scaled, 
                           labels_row=lab,
                           cluster_rows=TRUE, 
                           show_rownames=show_rownames,
                           #scale='row', 
                           cluster_cols=cluster_cols,
                           
                           annotation_col=df_ord, 
                           
                           clustering_method = 'ward.D2')
  
  

  #dev.off()
  
  return(my_pheatmap)
}


