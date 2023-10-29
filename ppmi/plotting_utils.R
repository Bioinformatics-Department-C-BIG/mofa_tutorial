


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
  dim(df_ord)
  show_rownames=TRUE
  jpeg(fname, width=10*100, height=7*100, res=100)
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
  dev.off()
  #ggsave(fname, width=7, height=7, dpi=300)
  
  ComplexHeatmap::pheatmap(hm_scaled, 
                           labels_row=lab,
                           cluster_rows=TRUE, 
                           show_rownames=show_rownames,
                           #scale='row', 
                           cluster_cols=cluster_cols,
                           
                           annotation_col=df_ord, 
                           
                           clustering_method = 'ward.D2')
  
  
  # ggsave(fname, width=10, height=7)
  


  return(my_pheatmap)
}





#remove_cn=FALSE
#order_by_hm='COHORT'

#cluster_cols=TRUE


#groups_kmeans3$cluster
#sel_samples=names(which(groups_kmeans3$cluster==3))

#sel_samples
#mt<-colData(vsd_filt)
#table(mt[mt$PATNO %in% sel_samples, 'NHY'])
#order_by_hm=c('PATNO_EVENT_ID')



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
  #vsd_filt=vsd
  
 
  sigGenes = make.names(sigGenes)
  rownames(vsd_filt) = make.names(rownames(vsd_filt))
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
  
  ## turn to sumbol? 
  
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
  
  
  dim(vsd_filt_genes)
  
  ### Add the annotations 
  meta_single<-colData(vsd_filt_genes)
  
  
  
  #get_symbols_vector
  ## HEATMAP OPTIONS 
  cluster_cols=cluster_cols;   
  
  #colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
  fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
                filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,remove_cn ,'.jpeg')
  
  hm<-assay(vsd_filt_genes); dim(hm);dim(df)
  
  df_ord=df
  hm_ord=hm
  
  ### SCALE!! 
  hm_scaled <- as.matrix(hm_ord) %>% t() %>% scale() %>% t()
  #jpeg(fname, width=2000, height=1500, res=200)
  library(ggplot2)
  if(process_mirnas){
    lab=rownames(rowData(vsd_filt_genes)) }else{
      lab=as.character(rowData(vsd_filt_genes)$SYMBOL)}

  
  
  filter_heatmap<-function(hm_scaled, df_ord,sel_pats=NULL, sel_pats_event_ids=NULL){
    # some inds are NA deal with it by subseting by patient 
    if (!is.null(sel_pats)){
      sel_pats=sel_pats[!is.na(sel_pats)]
      d_ind<-df_ord$PATNO %in% sel_pats
      
    }else{
      sel_pats_event_ids=sel_pats_event_ids[!is.na(sel_pats_event_ids)]
      d_ind<-df_ord$PATNO_EVENT_ID %in% sel_pats_event_ids
    }
    
    hm_scaled_filt<-hm_scaled[,d_ind]
    df_ord_filt<-df_ord[d_ind, ]
    return(list(hm_scaled_filt, df_ord_filt))
  }
  
  
  sel_pats<-df_ord[(df_ord$COHORT ==1), 'PATNO']; 
  res<-filter_heatmap(hm_scaled, df_ord, sel_pats)
  hm_scaled=res[[1]]; df_ord=res[[2]];
  dim(hm_scaled); dim(df_ord)
  
  
  show_rownames=TRUE
  #jpeg(fname, width=10*100, height=7*100, res=300)
  
  
  
  
  ######### ORDER OF COLUMNS/SAMPLES ####
  sel_pats_event_ids<-df_ord[(df_ord$EVENT_ID == 'BL'), 'PATNO_EVENT_ID']; 
  resV8<-filter_heatmap(hm_scaled, df_ord, sel_pats=NULL,sel_pats_event_ids = sel_pats_event_ids )
  hm_scaled_BL=resV8[[1]]; df_ord_BL<-resV8[[2]]
  
  
  
  sel_pats_event_ids<-df_ord[(df_ord$EVENT_ID == 'V08'), 'PATNO_EVENT_ID']; 
  resV8<-filter_heatmap(hm_scaled, df_ord, sel_pats=NULL,sel_pats_event_ids = sel_pats_event_ids )
  hm_scaled_v08=resV8[[1]]; df_ord_V08<-resV8[[2]]
  

  dend_v08<-as.dendrogram(hclust(dist(t(hm_scaled_v08))));
  order_V08_labels<-labels(dend_v08)
  order_v08<-seq(1:length(order_V08_labels)); names(order_v08)<-gsub('\\_.*','', order_V08_labels)
  
  df_ord$pat_order<-order_v08[match(df_ord$PATNO,names(order_v08)  )]
  
  
  new_ord<-order(df_ord[,'EVENT_ID'], df_ord[,'pat_order'])
  df_ord<-df_ord[new_ord ,]
  
  df_ord$COHORT
  dim(hm); dim(df_ord)
  hm_scaled<-hm_scaled[,new_ord]
  
  ####
  
  
  ## Put a blank row in each one 
  split_time<-as.numeric(cumsum(table(df_ord$EVENT_ID)))
  
  df_ord<-df_ord[,!colnames(df_ord) %in% c('PATNO_EVENT_ID', 'PATNO', 'pat_order')]
  
  if (!process_mirnas){
    lab<-get_symbols_vector(rownames(hm_scaled))
    
  }else{
    lab=rownames(hm_scaled)
  }
  # symbs
  
  graphics.off()
  
  na_cols<-apply(df_ord, 2, function(x){nas<-is.na(x);
  all(nas)})
  df_ord=df_ord[, !na_cols]
  
  dim(hm_scaled)
  dim(df_ord)
  
 
  
  df_ord[(df_ord$COHORT ==1),]
  rownames(df_ord)<-colnames(hm_scaled)
  
  my_pheatmap<-pheatmap(hm_scaled, 
                        labels_row=lab,
                        cluster_rows=TRUE, 
                        show_rownames=show_rownames,
                        gaps_col =  split_time,
                        cluster_cols=FALSE,
                        
                        annotation_col=df_ord, 
                        
                        clustering_method = 'ward.D2'
  )
  
  #dev.off()
  my_pheatmap
  
  
  
  # 1. SCALE
  # 2. Split time
  # 3. Cluster genes within groups? 
  library(ComplexHeatmap)
  library(circlize)
  
  # if only 1 timepoint
  
  hm_scaled_col = clip_outliers_times(hm_scaled_v08, x_times = 1.7)
  #hm_scaled_col = clip_outliers(hm_scaled_v08)
  
 # hm_scaled_col=hm_scaled
  f1 = colorRamp2(seq(min(hm_scaled_col), max(hm_scaled_col), length = 3), c("blue", "#EEEEEE", "red"))
  

  #factor_labels<-combined_sig_genes[match(rownames(hm1), combined_sig_genes$symbol), 'id']
  
  
  ### NOW SPLIT TO TWO heatmaps 
  
  df1<-df_ord_V08; hm1<-hm_scaled_v08
  hm_timepoint<-function(df1, hm1){
    #'
    #'
    #'
    df1<-df1[, !colnames(df1) %in% c('PATNO_EVENT_ID', 'PATNO')]
    ha1<-ComplexHeatmap::HeatmapAnnotation(df=df1)
    hm_t<-ComplexHeatmap::Heatmap(hm1, 
                                    top_annotation = ha1, 
                                    col=f1)
    return(hm_t)
  }
  
  ht_V08<-hm_timepoint(df_ord_V08, hm_scaled_v08)
  ht_BL<-hm_timepoint(df_ord_BL, hm_scaled_BL)
  
  #if (!is.null(df_ord_V06)){
  #  ht_V06<-hm_timepoint(df_ord_V06, hm_scaled_V06)
    
    
  #}
  
  # ggsave(fname, width=7, height=7, dpi=300)
  anno<-rowAnnotation(labels = anno_text(lab, which = "row"))

  
  
  
  hboth<-ht_V08+ht_BL+anno
  hboth<-ht_V08+anno
  
  factor_labels=NULL
  if (!is.null(factor_labels)){
    hm_draw<-draw(hboth, row_km = 1, row_split =factor_labels, cluster_rows = TRUE)
    
  }else{
    hm_draw<-draw(hboth, row_km = 1,cluster_rows = TRUE)
    
  }
  
  
  #dev.off()
  

  
  return(hm_draw)
}


#BiocManager::install('ComplexHeatmap')
