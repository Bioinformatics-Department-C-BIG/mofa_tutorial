

source(paste0('ppmi/setup_os.R'))
library(org.Hs.eg.db)

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))


### disconnect from mofa and other scripts 
VISIT=c('BL','V04', 'V06',  'V08');
VISIT=c('BL', 'V06' ,'V08');
VISIT=c('BL', 'V08');
VISIT=c('BL','V04', 'V06',  'V08');

process_mirnas=FALSE


source(paste0(script_dir,'ppmi/deseq_analysis_setup.R')) # this will load the required vsd 
source(paste0(script_dir,'ppmi/plotting_utils.R'))
source(paste0(script_dir,'ppmi/time_utils.R'))



write.csv(deseq2Results, paste0(outdir_s, '/results.csv'))

deseq2ResDF <- as.data.frame(deseq2Results)
deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
deseq2ResDF$abslog2pval<-abs(deseq2ResDF$log2pval)
min(deseq2ResDF$padj)
write.csv(deseq2ResDF, paste0(outdir_s, '/results_df.csv'))




if (process_mirnas){
  sel_view='miRNA'
  
}else{
  sel_view='RNA'
}
### Up to here output can be used for other reasons
##


RUN_DESEQ_ANALYSIS=FALSE





library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# Coerce to a data frame
#BiocManager::install('org.Hs.eg.db')
###  TODO: turn to gene symbols 
if (!process_mirnas){
  
  ens <- rownames(deseq2ResDF)
  ens<-gsub('\\..*', '',ens)
  symbols_ordered<-get_symbols_vector(ens)
  deseq2ResDF$SYMBOL<-symbols_ordered
  rowData(vsd)$SYMBOL=symbols_ordered
  rownames(deseq2ResDF)<-ens
  rownames(vsd)<-ens
  
  
}

log2fol_T<-0.25
padj_T<-.0005

deseq2ResDF_strict<-mark_significant(deseq2ResDF, padj_T, log2fol_T)
deseq2ResDF_strict<-mark_significant(deseq2ResDF, padj_T, log2fol_T)

signif_genes_strict<-rownames(deseq2ResDF_strict[!is.na(deseq2ResDF_strict$significant),])
signif_genes_strict

####### MOFA deseq2  

padj_T_hv<-0.05
log2fol_T_hv<-0.1

### this is also done later on -- save from there? 
deseq2ResDF$mofa_sign<- ifelse(deseq2ResDF$padj <padj_T_hv & abs(deseq2ResDF$log2FoldChange) >log2fol_T_hv , "Significant", NA)


#deseq2ResDF$mofa_sign

signif_genes<-rownames(deseq2ResDF[!is.na(deseq2ResDF$mofa_sign),])
length(signif_genes)

vsd_mat <- assay(vsd)
log2fol_T_overall<-0.1
padj_T_overall<-.05
deseq2ResDF
deseq2ResDF<-mark_significant(deseq2ResDF, padj_T_overall, log2fol_T_overall, outdir_single = outdir_s)





######## HEATMAP #######
  vsd_filt=vsd
 
  #Convert the VST counts to long format for ggplot2
  library(reshape2)
  library(pheatmap)
  library('pheatmap')
  
  
  df<-vsd_filt$COHORT
  assay(vsd_filt)
  
  
  colDataToPlot<-c('NP1RTOT','NP2_TOT', 'rigidity', 'td_pigd_old_on', 'moca' , 'RBD_TOT', 'NP3_TOT')
  colDataToPlot<-c('NP2_TOT', 'td_pigd_old_on',  'RBD_TOT', 'NP3_TOT', 'con_putamen', 'con_putamen_V10')
  colDataToPlot<-c('updrs3_score', 'td_pigd_old_on',  'RBD_TOT', 'updrs2_score', 'con_putamen', 
                   'updrs3_score_diff_V10', 'moca_diff_V10',  'updrs2_score_diff_V10',  'updrs3_score_V12', 'NP3_TOT_V16', 'con_putamen_V10', 
                   'con_putamen_diff_V10', 'PUTAMEN_R_V10', 
                   'mean_striatum_V10', 'abeta', 
                   'hi_putamen_V10')
  
  colDataToPlot<-c('updrs3_score','updrs3_score_V12', 'NP3_TOT_V16',   'updrs2_score', 
                   'con_putamen','con_putamen_V10', 'NP3_TOT_diff_V16', 'NP2_TOT_diff_V16',
                   'abeta', 'RBD_TOT_diff_V16', 
                   'hi_putamen_diff_V10', 'scopa')
  
  

  
  # TODO:  choose off 
  # DEFINE 
  # 1. SETTINGS 
  # 2. Clusters 
  # 3. factors and features 
  
  
  # 1. Annotation 
  
  y='NP2PTOT_LOG'
  y='NP3TOT_LOG'
    y='NP3TOT_LOG'
  y='sft'

  df_all<-fetch_metadata_by_patient_visit(vsd$PATNO_EVENT_ID , combined=combined_bl_log)
  df_all<-fetch_metadata_by_patient_visit(vsd$PATNO_EVENT_ID , combined=combined_bl_log)
  diff_variables= c('updrs2_score_diff_V12', 'NP2PTOT_diff_V14')
  colData_change<-c('updrs3_score', 'con_putamen', 'hi_putamen', 'updrs2_score', 'moca')
  sm<-MOFAobjectPD_sel@samples_metadata
  sm$NP3TOT_LOG_clust
  df_all<-cbind(df_all,sm[match(df_all$PATNO_EVENT_ID, sm$PATNO_EVENT_ID ),c(diff_variables, 'INEXPAGE') ])
  y='NP3TOT_LOG'
  y_clust<-paste0(y, '_clust')
  #y_diff=paste0(y, '_diff_V12')
  
  choose<-c(y,diff_variables,'scopa', "COHORT", "SEX", 'AGE', 'NHY','PATNO', 'EVENT_ID','PATNO_EVENT_ID','PDMEDYN')
  df<-as.data.frame(df_all[,choose])


  
  # 2.  CHOOSE FACTORS 
  # 2. and fetch features either from time diff or other diff
  heatmap_factors=which(all_fs_diff[,y]);heatmap_factors
  

  ### Plot MOFA too
 # if (length(VISIT)>1){
    remove_cn=FALSE
    cluster_cols=TRUE
    cluster_id = c(1,2)
    #sel_samples=names(which(groups_kmeans3$cluster%in% c(1,2)))

    
    #sel_samples
    mt<-colData(vsd)
    order_by_hm=c('PATNO_EVENT_ID')
    heatmap_factors
    
    ### TODO: source from molecular trajectories 
    #sigGenes_wil<-all_sig_genes2[grepl('ENS', all_sig_genes2$symbol),]
    #sigGenes=sigGenes_wil$symbol
    
    
    # 2. FACTORS   
    #sel_view
    # THE LOGIC IS to get the top genes from each factor 
    # TODO: what happens if one gene is in the top of two factors? get the highest rank/ normalized weight
    if (sel_view=='RNA'){top_fr=0.02}else{top_fr=0.05}
    top_fr
    f_features=concatenate_top_features(MOFAobject=MOFAobject,view=sel_view, factors_all=heatmap_factors,top_fr=top_fr)

    f_features$feature<-gsub('\\..*', '',f_features$feature)
    f_features=f_features[!duplicated(f_features$feature),]
    f_features
    
    factor_labels<-f_features$Factor; names(factor_labels)<-f_features$feature
    table(factor_labels)
    
    sigGenes_f=f_features$feature
    rnas_time<-get_de_rnas_union(sel_vis)
  rnas_time_all = unique(unlist(rnas_time))
  rnas_time_all = gsub('\\..*','', rnas_time_all)
  rnas_time_all
    sigGenes2<-intersect(sigGenes_f, rnas_time_all)
    sigGenes2

  plot_heatmap_time_clusters(vsd_filt=vsd,  sigGenes = sigGenes2 ,  df=df, remove_cn=FALSE )


    #### Which clusters? add to config
    clust_name<-paste0(y)
    clust_for_metric<-all_clusts_mofa[[clust_name]]
 


    df$cluster_m<-as.factor(sm$NP3TOT_LOG_clust[match(df$PATNO, sm$PATNO)] )
    vsd$cluster_m<-as.factor(sm$NP3TOT_LOG_clust[match(vsd$PATNO, sm$PATNO)] )
    vsd$cluster_m
    graphics.off()

    colnames(df)
    colDataToPlot<-colDataToPlot[colDataToPlot%in% colnames(df)]
    df[, c(colDataToPlot)]<-sapply(df[, colDataToPlot], as.numeric)
    df
    #df[df$EVENT_ID=='V08', 'con_putamen']=df[df$EVENT_ID=='V08', 'con_putamen_V10']
    df$AGE<-as.numeric(df$AGE)

    
    # todo: add the change as well
    
    ## TODO: symbols vector 
    
 

    
    # filter out samples here
    ## PROBLEM: cannot filter 

    sm_pd<-MOFAobjectPD_sel@samples_metadata
    
    sm_pd$NP3TOT_clust
    sel_clusts<-c(1,3)
    sm_pd[sm_pd[,clust_name] %in% sel_clusts,]
    sel_samples=sm_pd[sm_pd[,paste0(y, '_clust')] %in% sel_clusts,]$PATNO
#sel_samples

    draw_all_times=TRUE; wf<-100
    sm_pd=MOFAobjectPD@samples_metadata
    sel_cluster_ids=c(1,2,3)
    sel_cluster_ids=c(1,2)
    
    sel_cluster_ids_s=paste0(sel_cluster_ids, collapse='_')
    
    fname=paste0(outdir_s, '/heatmap_time', '_f_',paste(heatmap_factors, collapse='_'),'clust_', sel_cluster_ids_s,draw_all_times, '.jpeg')
    #  fname=paste0(outdir_s, '/heatmap/heatmap_time', '_f_',factor, '.jpeg')
    
    graphics.off()
    df$cluster_m
    my_pheatmap<-plot_heatmap_time(vsd_filt=vsd, sigGenes = sigGenes  ,  df=df, remove_cn=FALSE,
                                   show_rownames = show_rownames,cluster_cols = TRUE, sel_samples=sel_samples, 
                                   factor_labels=factor_labels,draw_all_times = draw_all_times)
    
    

    jpeg(fname, width=30*wf, height=12*200, res=200)
    
    my_pheatmap
    dev.off()
      
      


  
  



    
    
    #### DESEQ stages?? 
    
    # 1. match the clusters 
    # 2. deseq the clusters
    colData(se_filt)










