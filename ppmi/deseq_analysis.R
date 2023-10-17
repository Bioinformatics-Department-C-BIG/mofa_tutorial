

source(paste0('ppmi/setup_os.R'))

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))


### disconnect from mofa and other scripts 
#VISIT=c('BL','V04', 'V06',  'V08');
process_mirnas=TRUE

source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))
source(paste0(script_dir,'ppmi/plotting_utils.R'))



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

VISIT
#symbols[dup_ind]

outdir_s

deseq2ResDF$SYMBOL

deseq2ResDF


log2fol_T<-0.25
padj_T<-.005

deseq2ResDF_strict<-mark_significant(deseq2ResDF, padj_T, log2fol_T)
deseq2ResDF_strict<-mark_significant(deseq2ResDF, padj_T, log2fol_T)

min(deseq2ResDF$padj)

####### MOFA deseq2  

padj_T_hv<-0.05
log2fol_T_hv<-0.1

### this is also done later on -- save from there? 
deseq2ResDF$mofa_sign<- ifelse(deseq2ResDF$padj <padj_T_hv & abs(deseq2ResDF$log2FoldChange) >log2fol_T_hv , "Significant", NA)

#deseq2ResDF$mofa_sign

signif_genes<-rownames(deseq2ResDF[!is.na(deseq2ResDF$mofa_sign),])
length(signif_genes)

vsd_mat <- assay(vsd)

###TODO: Move this to the mofa file 


### TODO: ADD SIGNIFICANCE thresholds in the output file!! 
#for (most_var in c(0.05, 0.5)){
#  for (most_var in c(0.05, 0.1,0.15,0.2,0.25,0.3,  0.9,0.75,0.5)){
    for (most_var in c(0.05,0.1, 0.5, 0.9, 0.2,0.3, 0.35, 0.4, 0.45, 0.75)){
    

  param_str_tmp<-paste0(prefix, VISIT_S, '_',most_var ,'_', min.count, '_coh_', sel_coh_s, '_', sel_subcoh_s )
  highly_variable_outfile<-paste0(output_files, param_str_tmp,'_highly_variable_genes_mofa.csv')
  highly_variable_sign_outfile<-paste0(output_files, param_str_tmp,'_highly_variable_genes_mofa_signif.csv')
  # TODO: %features should stay the same after filter 
  highly_variable_genes_mofa<-selectMostVariable(vsd_mat, most_var)
  highly_variable_sign_genes_mofa<-highly_variable_genes_mofa[rownames(highly_variable_genes_mofa) %in%  signif_genes,]
  
  
  write.csv(highly_variable_genes_mofa, highly_variable_outfile);
  
  write.csv(highly_variable_sign_genes_mofa, highly_variable_sign_outfile)

  
  highly_variable_sign_outfile
  rownames(highly_variable_genes_mofa)
  
  highly_variable_outfile
  
  
} 
dim(highly_variable_sign_genes_mofa)



log2fol_T_overall<-0.1
padj_T_overall<-.05
deseq2ResDF
deseq2ResDF<-mark_significant(deseq2ResDF, padj_T_overall, log2fol_T_overall, outdir_single = outdir_s)


num_de_genes<-length(which(!is.na(deseq2ResDF$significant)))

#hist(-log10(deseq2ResDF$padj))
#hist(deseq2ResDF$log2FoldChange)
#hist(deseq2ResDF$log2pval[abs(deseq2ResDF$log2pval)>0.15])
#
#
#
#min(deseq2ResDF$log2pval, na.rm = TRUE)





#### PCA plots
### make the plot with higly variable AND significant genes 




pca_data_list<-list(t(highly_variable_sign_genes_mofa),t(highly_variable_genes_mofa) )

pc_ind_ps<-list()
for (i in c(1,2)){
      pca_data<-pca_data_list[[i]]
  

      pca_pars<-paste0('_signif_', i)
      
      dim(highly_variable_sign_genes_mofa)
      meta_d<-colData(se_filt)
      coh_d<-meta_d$COHORT_DEFINITION
      
      pca.data <- PCA(pca_data,
                      scale.unit = TRUE, graph = FALSE
                      )
      
      fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))
      
      
      pc_ind_p<-fviz_pca_ind(pca.data,col.ind = meta_d$COHORT_DEFINITION, 
                   label='none')
      pc_ind_ps[[i]]=pc_ind_p
      ggsave(paste0(pca_files, 'individuals', pca_pars, '.jpeg'), width=6, height = 6)
      
      ### graph of variables 
      
      
      
      fviz_pca_var(pca.data,
        col.var="contrib")+
        scale_color_gradient2(low="blue", mid="white", 
                              high="red", midpoint=1.5)+theme_bw()
      
      
      #### TODO: calculate PCs for ALL the datasets 
      
      
      ggsave(paste0(pca_files, 'variables', pca_pars,'.jpeg'), width=6, height = 6)
}
#### Make plots 
# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))








###### PLOTS of DE genes 

## Plot the results similar to DEseq2

# Let's add some more detail

limits<-as.numeric(max(abs(deseq2ResDF$log2FoldChange)))
limits
p_log_plot<-ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-2, 2), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, 
                               linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() +
  geom_density_2d(colour="black", size=0.5)+ 
  geom_hline(yintercept=c(0.15, -0.15),linetype="dashed", color = "red")+
  ylim(-(limits+0.1),(limits+0.1))
  

p_log_plot

logplot_f<-paste0(outdir_s, '/logfold_padjusted_plot.jpeg')
ggsave(logplot_f,width=8, height=6 )


#### NORMALIZED COUNTS
top_gene<-rownames(deseq2ResDF[which.min(deseq2ResDF$padj),])
#top_gene

contrasts<-c('Condition', 'EVENT_ID')
dds$Condition<-dds$COHORT
  
# Extract counts for the gene otop2

##### 
# Compute normalization factors and vst 



#dds <- estimateSizeFactors(dds,)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function
#vsd <- varianceStabilizingTransformation(dds)
graphics.off()
run_heatmap=TRUE


### TODO: MAKE A FUNCTION 
if (run_heatmap){
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
  order_by_metric<-'abslog2pval'
  order_by_metric<-'log2pval'
  # bring the low padj to the top!!
  order_by_metric<-'abslog2pval'
  order_by_metric<-'abslog2'
  order_by_metric<-'abslog2FC'
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
  oSigGenes
  
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
  
  
  
  
  
 # detach('ComplexHeatmap',unload=TRUE)
         
  
   
  
  df<-vsd_filt$COHORT
  assay(vsd_filt)
  
  
  colDataToPlot<-c('NP1RTOT','NP2_TOT', 'rigidity', 'td_pigd_old_on', 'moca' , 'RBD_TOT', 'NP3_TOT')
  colDataToPlot<-c('NP2_TOT', 'td_pigd_old_on',  'RBD_TOT', 'NP3_TOT')
  
  df<-as.data.frame(colData(vsd_filt)[,c( "SEX", 'AGE', 'NHY','PATNO', 'EVENT_ID','PDMEDYN', colDataToPlot, "COHORT")])
  # if clusters exist 
  
  
  # clusters_single
#  df$cluster_s<-factor(clusters_single$cluster[match(colData(vsd_filt)$PATNO_EVENT_ID, names(clusters_single$cluster ))])
  #### Add different clustering? 
 # df$cluster_m<-factor(clusters_mofa$cluster[match(colData(vsd_filt)$PATNO_EVENT_ID, names(clusters_mofa$cluster ))])
  #df$cluster_m_34<-factor(clusters_mofa_34$cluster[match(colData(vsd_filt)$PATNO_EVENT_ID, names(clusters_mofa_34$cluster ))])
  
  
  
  ws_top_bottom=select_top_bottom_perc(view='RNA', factors=c(3,4))

  graphics.off()
  
  

    my_pheatmap<-plot_heatmap(vsd_filt=vsd_filt, sigGenes = ws_top_bottom  ,  df=df, remove_cn=FALSE,
                            show_rownames = show_rownames,cluster_cols = TRUE, sel_samples = NULL )
  
    
    remove_cn=FALSE
    order_by_hm='COHORT'
    
    cluster_cols=TRUE
    
    
    groups_kmeans3$cluster
    sel_samples=names(which(groups_kmeans3$cluster==3))
    
    #sel_samples
    mt<-colData(vsd_filt)
    table(mt[mt$PATNO %in% sel_samples, 'NHY'])
    order_by_hm=c('PATNO_EVENT_ID')
    sigGenes<-most_sig_over_time$symbol[1:40]
    graphics.off()
  ### Plot MOFA too
  if (length(VISIT)>1){
    
    ## TODO: symbols vector 
    my_pheatmap<-plot_heatmap_time(vsd_filt=vsd_filt, sigGenes = sigGenes  ,  df=df, remove_cn=FALSE,
                                   show_rownames = show_rownames,cluster_cols = TRUE, sel_samples=sel_samples )
    
    my_pheatmap
  }
    

  


  
  
  #plot(hclust(dists))
}



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









