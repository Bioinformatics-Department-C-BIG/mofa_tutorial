

source(paste0('ppmi/setup_os.R'))

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))


### disconnect from mofa and other scripts 
VISIT=c('BL','V04', 'V06',  'V08');
VISIT=c('BL', 'V06' ,'V08');
VISIT=c('BL','V04', 'V06',  'V08');
VISIT=c('BL', 'V08');

#process_mirnas=FALSE

source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))
source(paste0(script_dir,'ppmi/plotting_utils.R'))
source(paste0(script_dir,'ppmi/time_utils.R'))



write.csv(deseq2Results, paste0(outdir_s, '/results.csv'))

deseq2ResDF <- as.data.frame(deseq2Results)
deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
deseq2ResDF$abslog2pval<-abs(deseq2ResDF$log2pval)
min(deseq2ResDF$padj)
write.csv(deseq2ResDF, paste0(outdir_s, '/results_df.csv'))


### Up to here output can be used for other reasons
##

ddsSE@design

deseq2Results@metadata

RUN_DESEQ_ANALYSIS=FALSE

# 
run_ma<-FALSE
if (run_ma){
  jpeg(paste0(outdir_s, '/MA_plot_results.jpeg'))
  plotMA(deseq2Results)
  dev.off()
}



#, ylim=c(-1,10), xlim=c(0,5))




library(ggplot2)


library(scales) # needed for oob parameter
library(viridis)
#install.packages('viridis')
#install.packages('gridExtra')


# Coerce to a data frame

### TODO: what is lfcShrink? 

library(org.Hs.eg.db)
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
write.csv(deseq2ResDF, paste0(outdir_s, '/results_df.csv'))

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
  ### INPUT FOR HEATMAP IS VSD
  deseq2VST <- as.data.frame( assay(vsd_filt))
  rownames(deseq2VST)<-rownames(vsd)
  deseq2VST$Gene <- rownames(deseq2VST)

  
  #deseq2ResDF$padj 
  # Keep only the significantly differentiated genes where the fold-change was at least 3
  log2fol_T_hm<-0.1
  padj_T_hm<-.05
  deseq2ResDF$abslog2FC<-abs(deseq2ResDF$log2FoldChange)
  sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= padj_T_hm & abs(deseq2ResDF$log2FoldChange) > log2fol_T_hm,])
  deseq2ResDF$Gene<-rownames(deseq2ResDF)

  order_by_metric<-'padj_reverse'
  
  quant=0.9
  mean_expr_T<-quantile(deseq2ResDF$baseMean, quant)
  # DEFINE THE reverse of padj
  deseq2ResDF$padj_reverse<--deseq2ResDF$padj
  deseq2ResDF$padj_reverse

  ### Filter the significant 
  oSigGenes<-deseq2ResDF[deseq2ResDF$padj <= padj_T_hm  & abs(deseq2ResDF$log2FoldChange) > log2fol_T_hm, ] ; dim(oSigGenes)
  if (dim(oSigGenes)[1]<1){
    log2fol_T_hm<-0
    padj_T_hm<-.05
    oSigGenes<-deseq2ResDF[deseq2ResDF$pvalue <= padj_T_hm  & abs(deseq2ResDF$log2FoldChange), ] ; dim(oSigGenes)
    
  }

  length(rownames(highly_variable_sign_genes_mofa))
  #View(deseq2ResDF)
  filter_highly_var=FALSE
  if (filter_highly_var){
    oSigGenes<-deseq2ResDF[rownames(deseq2ResDF) %in% rownames(highly_variable_sign_genes_mofa),]
    dim(oSigGenes)
    most_var_t=most_var
  }else{
    most_var_t=FALSE
  }
  ### also filter by mean > 0.75
  ## DO NOT FILTER BY MEAN 
  #oSigGenes<-oSigGenes[oSigGenes$baseMean>mean_expr_T,];dim(oSigGenes)
  
  ### Order using the order_by_metric 
  orderedSigGenes<-oSigGenes[order(-oSigGenes[,order_by_metric]),]
  orderedSigGenes
  n_sig_f='all'
  n_sig_f=1000; show_rownames=FALSE
  
  if (n_sig_f=='all'){
    n_sig=dim(orderedSigGenes)[1]
    
  }else{
    n_sig=n_sig_f
  }
  
  
  ## filter the top 50 significant genes 
  sigGenes <- orderedSigGenes$Gene[1:n_sig]
  deseq2VST_sign <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
  
  dim(deseq2VST)
  #Convert the VST counts to long format for ggplot2
  library(reshape2)
  library(pheatmap)
  
  deseq2VST_wide <- deseq2VST_sign
  deseq2VST_long <- melt(deseq2VST_wide, id.vars=c("Gene"))
  
  # Make a heatmap ## using the vst
  ## TODO add annotation
  

  
  
    
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
                   'hi_putamen_diff_V10')
  
  

  
  # TODO:  choose off 
  # DEFINE 
  # 1. SETTINGS 
  # 2. Clusters 
  # 3. factors and features 
  
  
  # 1. Annotation 
  df_all<-fetch_metadata_by_patient_visit(vsd$PATNO_EVENT_ID , combined=combined_bl_log)
  colData_change<-c('updrs3_score', 'con_putamen', 'hi_putamen', 'updrs2_score', 'moca')
  sm<-MOFAobject@samples_metadata
  df_all<-cbind(df_all,sm[match(df_all$PATNO_EVENT_ID, sm$PATNO_EVENT_ID ),c(diff_variables, 'INEXPAGE') ])
  y='NP2PTOT'
  y_clust<-paste0(y, '_clust')
  
  choose<-c( "SEX", 'AGE', 'NHY','PATNO', 'EVENT_ID','PATNO_EVENT_ID','PDMEDYN',y, "COHORT")
  df<-as.data.frame(df_all[,choose])


  
  # 2. 
  heatmap_factors=which(all_fs_diff[,y])
  

  ### Plot MOFA too
  if (length(VISIT)>1){
    remove_cn=FALSE
    cluster_cols=TRUE
    cluster_id = c(1,2)
    #sel_samples=names(which(groups_kmeans3$cluster%in% c(1,2)))

    
    #sel_samples
    mt<-colData(vsd)
    table(mt[mt$PATNO %in% sel_samples, 'NHY'])
    order_by_hm=c('PATNO_EVENT_ID')

    
    ### TODO: source from molecular trajectories 
    #sigGenes_wil<-all_sig_genes2[grepl('ENS', all_sig_genes2$symbol),]
    #sigGenes=sigGenes_wil$symbol
    
    
    # 2. FACTORS   
    sigGenes= combined_sig_genes_strict$symbol
    sigGenes= signif_genes_strict
    sigGenes= combined_sig_genes_strict$symbol
    
    ws_top_bottom=select_top_bottom_perc(view='RNA', factors=heatmap_factors, top_fr=.009 )
   
    ### TODO: tomorrow get the most variable or the most 
    
    #factor_labels<-combined_sig_genes[match(rownames(hm1), combined_sig_genes$symbol), 'id']
    #factor_labels<-combined_sig_genes[match(rownames(hm1), combined_sig_genes$symbol), 'id']
    
    
    ws_top_bottom<-gsub('\\..*', '',ws_top_bottom)
    sigGenes=ws_top_bottom
    sigGenes

    #### Which clusters? add to config
    clust_name<-paste0(y)
    clust_for_metric<-all_clusts_mofa[[clust_name]]
    df$cluster_m<-as.factor(clust_for_metric[match(df$PATNO_EVENT_ID, names(clust_for_metric))] )
    df$cluster_m
    graphics.off()

    colnames(df)
    
    df[, c(colDataToPlot)]<-sapply(df[, colDataToPlot], as.numeric)
    df
    df[df$EVENT_ID=='V08', 'con_putamen']=df[df$EVENT_ID=='V08', 'con_putamen_V10']
    df$AGE<-as.numeric(df$AGE)

    
    # todo: add the change as well
    
    ## TODO: symbols vector 
    
    
    fname=paste0(outdir_s, '/heatmap_time', '_f_',factor, '.jpeg')
    graphics.off()

    
    # filter out samples here
    ## PROBLEM: cannot filter 
    sel_samples=sm$PATNO_EVENT_ID
    my_pheatmap<-plot_heatmap_time(vsd_filt=vsd, sigGenes = sigGenes  ,  df=df, remove_cn=FALSE,
                                   show_rownames = show_rownames,cluster_cols = TRUE, sel_samples=NULL )
    
    
      jpeg(fname, width=10*150, height=12*200, res=200)
    
      my_pheatmap
      dev.off()
      
      
    
    
    }


  
  
  #plot(hclust(dists))




### Add Volcano plots 
### Compare to their results 



deseq2ResDF$SYMBOL

library('EnhancedVolcano')
if(process_mirnas){lab=rownames(deseq2ResDF) }else{lab=deseq2ResDF$SYMBOL}

mfc<-max(abs(deseq2ResDF$log2FoldChange))
pmax<-max(-log10(deseq2ResDF$padj), na.rm = TRUE)
pmax
xlim = c(-mfc-0.2,mfc+0.2)
ylim = c(0,pmax+0.2)
ylim = c(0,pmax-0.5)

ns_full<-table(se_filt$COHORT_DEFINITION)
ns<-paste0(rownames(ns_full)[1],' ', ns_full[1], '\n' ,names(ns_full)[2], ' ', ns_full[2])
ns
pvol<-EnhancedVolcano(deseq2ResDF,
                lab = lab,
                pCutoff = 10e-2,
                FCcutoff = 0.1,
                x = 'log2FoldChange',
                y = 'padj',
                
                
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
                title=''
               )


pvol

fname
fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT,'.jpeg')
ggsave(fname,pvol, width=4.5,height=7, dpi=300)


Padj_T_paths=0.05
padj_paths<-Padj_T_paths
pvalueCutoff=1
if (!process_mirnas){
  source('ppmi/RNAseq enrichment.R')
  
}else{
  source('ppmi/miRNA_seq_enrichment.R')
  
}









