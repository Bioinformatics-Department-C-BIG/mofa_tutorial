


### Complex heatmap 

  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  suppressPackageStartupMessages(library(circlize))
   suppressPackageStartupMessages(library('EnhancedVolcano'))


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
  #fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
  #              filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,remove_cn ,'.jpeg')
  
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

get_limits<-function(x){
        xmax<-max(x, na.rm = TRUE)
        xmin<-min(x,na.rm = TRUE)
        max_all<-max(c(abs(xmax), abs(xmin)))
        xmax
        return(c(-max_all, max_all))
}


create_heatmap_proteins<-function(hm,se_filt, fname,coldata_to_plot=c()){
  #' 
  #' simple heatmap for proteins 
  #' 
    
    df<-as.data.frame(colData(se_filt)[coldata_to_plot]); rownames(df)<-df$PATNO_EVENT_ID
    df$PATNO_EVENT_ID<-NULL
#rownames(df)
#colnames(hm)
    jpeg(fname, res = 200, width=3+log2(dim(hm)[2]), height=3+log2(dim(hm)[1]), units='in')

        my_pheatmap<-pheatmap(hm, 
                              #labels_row=lab,
                              cluster_rows=TRUE, 
                              show_rownames=TRUE,
                              scale='row', 
                              cluster_cols=cluster_cols,
                              annotation_col=df, 
                              clustering_method = 'complete', 
                              clustering_distance_cols = 'euclidean'
        )
        dim(hm)

        show(my_pheatmap)
        dev.off()
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





plotVolcano<-function(deseq2ResDF, se_filt, title='', xlim=NULL, lab=NULL,x='log2FoldChange', y= 'padj', FCcutoff=0.1,pCutoff=10e-2){
  #'
  #'
  #' Take a sumarized experiment and deseq results 
  #' @param deseq2ResDF deseq results dataframe 
  #'  @param summarized experiment 
  #' @param x metric to plot on x eg. log2FoldChange
  #' @param y padj
  #'
  #'
  #'
  
  
  #if(process_mirnas){lab=rownames(deseq2ResDF) }else{lab=deseq2ResDF$GENE_SYMBOL}
    

  mfc<-max(abs(deseq2ResDF[,x]))
  pmax<-max(-log10(deseq2ResDF[,y]), na.rm = TRUE)
  
  if (is.null(xlim)){ # if not supplied create it 

    xlim = c(-mfc-0.2,mfc+0.2)

  }
  ylim = c(0,pmax+0.2)
  ylim = c(0,pmax-0.5)
  
  ns_full<-table(se_filt$COHORT_DEFINITION)
  ns<-paste0(rownames(ns_full)[1],' ', ns_full[1], '\n' ,names(ns_full)[2], ' ', ns_full[2])
  
  if (is.null(lab)){
    lab=rownames(deseq2ResDF)
  }

  pvol<-EnhancedVolcano(deseq2ResDF,
                        lab = lab,
                        pCutoff = pCutoff,
                        FCcutoff =FCcutoff,
                        x = x,
                        y = y,
                        
                        ## format 
                        pointSize = 2,
                        legendIconSize = 5,
                        labSize = 4,
                        
                        legendLabSize=16,
                        subtitleLabSize = 13,
                        axisLabSize=17,
                        colAlpha = 0.5,
                        
                        # legend positions 
                        # legendPosition = 'right',
                        
                        xlim=xlim, 
                        ylim=ylim, 
                        
                        subtitle=ns, 
                        title=title
  )
  
  
  pvol
  return(pvol)
}



plotVolcano_proteins<-function(results_de, se_filt, title='', xlim=NULL, lab=NULL){



}



plot_heatmap_time<-function(vsd_filt, sigGenes,  df,remove_cn=FALSE, show_rownames=TRUE, 
                            cluster_cols=FALSE, order_by_hm='COHORT', sel_samples, 
                            factor_labels=NULL, draw_all_times=FALSE){
  
  
  
  #'
  #' @param vsd_filt: annotation dataframe nsamples X ncoldata 
  #' @param hm: heatmap data nfeats X nsamples 
  #' @param sigGenes: select genes to plot description
  #' @param factor_labels if not NULL it splits the heatmap into rows of factors  description
  #' 
  #' 
  #' 
  #' 
 #ARRANGE
  #vsd_filt=vsd
  
  #@draw_all_times=TRUE
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
  #fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
  #              filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols,remove_cn ,'.jpeg')
  
  hm<-assay(vsd_filt_genes); dim(hm);dim(df)
  
  df_ord=df
  hm_ord=hm
  
  ### SCALE!! 
  hm_scaled <- as.matrix(hm_ord) %>% t() %>% scale() %>% t()
  #jpeg(fname, width=2000, height=1500, res=200)
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
  if (!is.null(sel_samples)){
    sel_pats<-sel_pats[sel_pats%in%sel_samples]
  }
  res<-filter_heatmap(hm_scaled, df_ord, sel_pats)
  hm_scaled=res[[1]]; df_ord=res[[2]];
  dim(hm_scaled); dim(df_ord)
  

  
  show_rownames=TRUE
  #jpeg(fname, width=10*100, height=7*100, res=300)
  
  
  
  
#### FILTER THE TWO TIME POINTS
  sel_pats_event_ids<-df_ord[(df_ord$EVENT_ID == 'BL'), 'PATNO_EVENT_ID']; 
  resV8<-filter_heatmap(hm_scaled, df_ord, sel_pats=NULL,sel_pats_event_ids = sel_pats_event_ids )
  hm_scaled_BL=resV8[[1]]; df_ord_BL<-resV8[[2]]
  
  
  
  sel_pats_event_ids<-df_ord[(df_ord$EVENT_ID == 'V08'), 'PATNO_EVENT_ID']; 
  resV8<-filter_heatmap(hm_scaled, df_ord, sel_pats=NULL,sel_pats_event_ids = sel_pats_event_ids )
  hm_scaled_v08=resV8[[1]]; df_ord_V08<-resV8[[2]]
  
  
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
  

  
  
  # 1. SCALE
  # 2. Split time
  # 3. Cluster genes within groups? 

  
  # if only 1 timepoint
  
  if (draw_all_times){
    hm_scaled_col = clip_outliers_times(hm_scaled, x_times = 1.7)
    
  }else{
    hm_scaled_col = clip_outliers_times(hm_scaled_v08, x_times = 1.7)
    
  }
  #hm_scaled_col = clip_outliers(hm_scaled_v08)
  
 # hm_scaled_col=hm_scaled
  f1 = colorRamp2(seq(min(hm_scaled_col), max(hm_scaled_col), length = 3), c("blue", "#EEEEEE", "red"))
  
  ### NOW SPLIT TO TWO heatmaps 
  
  df1<-df_ord_V08; hm1<-hm_scaled_v08
  remove_from_hm<-c('PATNO_EVENT_ID', 'PATNO', 'COHORT')
  hm_timepoint<-function(df1, hm1, tp_name, cluster_columns=TRUE){
    #'
    #'
    #'
    df1<-df1[, !colnames(df1) %in% remove_from_hm]
    ha1<-ComplexHeatmap::HeatmapAnnotation(df=df1)
    hm_t<-ComplexHeatmap::Heatmap(hm1, 
                                    top_annotation = ha1, 
                                    col=f1, 
                                  cluster_columns=cluster_columns,
                                  
                                  name=tp_name)
    return(hm_t)
  }
  
  ht_V08<-hm_timepoint(df_ord_V08, hm_scaled_v08, tp_name='V08')
  ### padd V08 order to BL

  order_v08<-colnames(ht_V08@matrix)[column_order(ht_V08)]
  order_v08<-gsub('\\_.*','', order_v08)
  
  df_ord_BL_match<-df_ord_BL[match(order_v08,df_ord_BL$PATNO ),]
  bl_names<-gsub('\\_.*','',colnames(hm_scaled_BL))
  hm_scaled_BL_match<-hm_scaled_BL[,match(order_v08,bl_names )]
  
  
  ht_BL<-hm_timepoint(df_ord_BL_match, hm_scaled_BL_match,tp_name='BL', cluster_columns=FALSE)
  
  #if (!is.null(df_ord_V06)){
  #  ht_V06<-hm_timepoint(df_ord_V06, hm_scaled_V06)
    
    
  #}
  
  # ggsave(fname, width=7, height=7, dpi=300)
  anno<-rowAnnotation(labels = anno_text(lab, which = "row"))

  
  
  
  
  if (draw_all_times){
    hboth<-ht_BL+ht_V08+anno
    
  }else{
    hboth<-ht_V08+anno
    
  }
  
  
  if (!is.null(factor_labels)){
    hm_draw<-draw(hboth, row_km = 1, row_split =factor_labels, cluster_rows = TRUE, main_heatmap='V08')
    
  }else{
    hm_draw<-draw(hboth, row_km = 1,cluster_rows = TRUE,  main_heatmap='ht_V08')
    
  }
  
  
  #dev.off()
  

  
  return(hm_draw)
}


#BiocManager::install('ComplexHeatmap')

median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}





plot_molecular_trajectories<-function(merged_melt_filt_most_sig){
  #####
  #####
  #''
  
  
  ggplot(data = merged_melt_filt_most_sig, aes(x = VISIT, y = value, fill=kmeans_grouping)) + 
    #geom_point(aes(col=VISIT), size = 2) +
    #geom_line(aes(group=PATNO),  col= 'grey') +
    # subgroup should be in the fill parameter!!! 
    geom_boxplot(aes(x=VISIT, fill=kmeans_grouping ))+
    scale_color_viridis_d(option='mako')+
    scale_fill_viridis_d(option='mako')+
    
    #geom_line(aes(group=patno), palette='jco') +
    #facet_wrap(. ~ symbol) +
    
    geom_signif(comparisons = list(c('BL', 'V08')),  
                map_signif_level=TRUE, 
                tip_length = 0, vjust=0.4)+
    
    facet_wrap(. ~ symbol, scales='free_y') +
    
    theme_bw() 
  ggsave(paste0(outdir, '/trajectories/boxplots_',factor,'_',  keep_all_feats,'_cl_fs_',factors_to_cluster_s, view,'_',group_cat,sel_cohort , '.jpeg'), 
         width=20, height=12)
  
 
  
}




boxplot_by_cluster_multiple<-function(met, clust_name, diff_variables_to_p, bn,  height=10/1.5, width=10,facet_rows=1, text=''){
  #'
  #' Create boxplot by cluster 
  #' 
  #' @param met
  #' @param bn: name of metric
  #'  
  #' 
  clust_metric<-gsub('_clust', '', clust_name)
  

  as.numeric(met[,clust_metric])
  met[,clust_metric ]<-as.numeric(met[,clust_metric])
  met<-met[!is.na(met[, clust_metric]),]
  #print(paste('Using subset of  ', dim(met)[1], ' patients'))
  freqs<-paste0('n=', paste0(table(met[, clust_name]), collapse = ', '))
  


  
  #### PROPORTIONS OF BINARY VARS
  tot_med<-as.matrix(table(met[,c(clust_name, "PDMEDYN")])); paste_med<-paste0('Med: ' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ))
  tot_med<-as.matrix(table(met[,c(clust_name, "SEX")])); paste_sex<-paste('SEX:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  tot_med<-as.matrix(table(met[,c(clust_name, "PDSTATE")])); paste_state<-paste('PD state:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  tot_med<-as.matrix(table(met[,c(clust_name, "td_pigd")])); paste_tdpigd<-paste('td/pigd:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  tot_med<-as.matrix(table(met[,c(clust_name, "NHY")])); paste_nhy<-paste('H&Y:' ,paste0(format(tot_med[,2]/ rowSums(tot_med), digits=2), collapse=',' ) )
  
  k_centers<-max(as.numeric(unique(met[!(met[, clust_name] %in% 'HC'), clust_name] )) , na.rm = TRUE)
  k_centers
  ## Add kruskal wallis to the total plot and separately to each one 
  ## 
  factors<-paste0(which(all_fs_diff[,clust_metric]), collapse=', ')


  met_diff<-met[,c( 'PATNO',clust_name,diff_variables_to_p)]

  met_diff_val=reshape::melt(met_diff, id.vars=c('PATNO', clust_name))
  met_diff_val[, clust_name] = as.factor(met_diff_val[, clust_name] )
  met_diff_val[, 'value'] = as.numeric(met_diff_val[, 'value'] )
#  num_log<-is.numeric(met_diff_val)

  #met_diff_val[,num_log]<-sapply(met_diff_val[,num_log], clipping_values)



    p<-ggplot(met_diff_val ,aes_string(x=clust_name , y='value'))+
    geom_boxplot(aes_string( x=clust_name,# color=clust_name, 
                             fill=clust_name, y='value'), alpha=0.9, 
                             outlier.shape = NA)+

  
    facet_wrap(~variable, scales='free',#, labeller=labeller(
      #variable=labeller_clinical
      nrow=facet_rows)+
    geom_pwc( tip.length = 0,
      method = "wilcox_test", label = "p.adj.signif", label.size = 2)+   
    scale_color_viridis_d(option="turbo")+
       scale_fill_viridis_d(option="turbo")+
  
    
    
  labs(#title = paste(y),  
         #subtitle=paste(freqs, '\n','Kruskal.wallis p.val', format(kw$p.value, digits=2)),
         caption = paste0( '\n', 'factors: ',factors, ',  ', freqs,',  ', 
                          paste_med, ',  ',
                          paste_sex, ' ',
                          #paste_state,  ',\n',
                          # paste_tdpigd,' '
                          text), 
         x=mt_kv[which(mt_kv[,1]==clust_metric),2])+
    guides(fill=guide_legend(title='PD subgroup' ), color=guide_legend(title='PD subgroup' ))+
    theme(text =element_text(size=15), axis.title.y=element_blank())+

    theme(strip.text = element_text(face = "bold"))
  #plot.title = element_text(size = 30, color = "green")
  p

  print(bn)
  ggsave(bn, dpi=200, width=width, height = height, units='in')
  graphics.off()
  ## TODO: WILCOX TEST BY GROUP
  
  
}

#x='month'
#add_patient_lines=FALSE

plot_molecular_trajectories_line<-function(merged_melt_filt_most_sig, x='month',add_patient_lines=FALSE, trajectory_fname ){
  #'
  #' @param merged_melt_filt_most_sig
  #' @param x axis
  #' @param add_patient_lines
  #'
  #'
 # nrow=4
  
  p<-ggplot(data = merged_melt_filt_most_sig, aes_string(x = x, y = 'value', 
                                                         fill='group', group='group', colour='group')) 
  if (add_patient_lines){
    p<- p+geom_line(aes_string(x = x, y = 'value', 
                               group='PATNO', colour='group' ),size=0.1, alpha=0.5)
  }
  
  
  p=p+ stat_summary(geom = "errorbar", fun.data = median_IQR, 
                    position=position_dodge(0), alpha=0.9, width=0.3)+
    # horizontal lines 
    stat_summary(fun = median, position=position_dodge(width=0), 
                 geom = "line", size = 0.9, alpha=0.9 ) + # , linetype='longdash' 
    scale_color_viridis_d(option='turbo')+
    facet_wrap(. ~ symbol, scales='free_y', 
               nrow = nrow) +
    
    #geom_signif(comparisons = list(c('BL', 'V08')), 
    #            map_signif_level=TRUE, 
    #            tip_length = 0, vjust=0.3)+
    
    labs(y='logCPM')+
    # legend(legend=c('Low', 'High'))+
    theme(strip.text = element_text(
      size = 12, color = "dark green", face="bold"), 
      axis.title.y =element_text(
        size = 12, color = "dark green", face="bold",), 
      axis.text.x = element_text(
        size = 12 ))+
    guides(fill=guide_legend(title='PD subgroup' ), color=guide_legend(title='PD subgroup' ))
  
  
  
  p
  #warnings()
  ggsave(trajectory_fname, 
         width=width, height=height, dpi = 300)
  
  graphics.off()
  
  
  return(p)
  
}



