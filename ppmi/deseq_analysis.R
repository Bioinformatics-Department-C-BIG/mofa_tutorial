

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




library(ggplot2)


library(scales) # needed for oob parameter
library(viridis)
#install.packages('viridis')
#install.packages('gridExtra')


# Coerce to a data frame

### TODO: what is lfcShrink? 

library(org.Hs.eg.db)
#BiocManager::install('org.Hs.eg.db')

# TODO: function!! 
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

###TODO: Move this to the mofa file 


### TODO: ADD SIGNIFICANCE thresholds in the output file!! 
#for (most_var in c(0.05, 0.5)){
#  for (most_var in c(0.05, 0.1,0.15,0.2,0.25,0.3,  0.9,0.75,0.5)){
#get_highly_variable_matrix(prefix='mirnas_', VISIT_S, MIN_COUNT_M, sel_coh_s , sel_subcoh_s , TOP_MN)


    for (most_var in c(0.05,0.1, 0.5, 0.9, 0.2,0.3, 0.35, 0.4, 0.45, 0.75)){
    

  param_str_tmp<-paste0(prefix, VISIT_S, '_',most_var ,'_', min.count, '_coh_', sel_coh_s, '_', sel_subcoh_s )
  highly_variable_outfile<-paste0(output_files, param_str_tmp,'_highly_variable_genes_mofa.csv')
  highly_variable_sign_outfile<-paste0(output_files, param_str_tmp,'_highly_variable_genes_mofa_signif.csv')
  # TODO: %features should stay the same after filter 
  highly_variable_genes_mofa<-selectMostVariable(vsd_mat, most_var)
  highly_variable_sign_genes_mofa<-highly_variable_genes_mofa[rownames(highly_variable_genes_mofa) %in%  signif_genes,]
  
  
  #write.csv(highly_variable_genes_mofa, highly_variable_outfile);
  
  #write.csv(highly_variable_sign_genes_mofa, highly_variable_sign_outfile)

  
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







deseq2ResDF$SYMBOL


# TODO: function take se filt and deseq 
plate

cluster_id = 1

se_filt=se_filt_all[[cluster_id]]
deseq2ResDF=deseq_all_groups[[cluster_id]]

pvol<-plotVolcano(deseq2ResDF, se_filt)


pvol
fname
fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT,'.jpeg')
fname<-paste0(outdir_s, '/EnhancedVolcano_edited_', prefix, VISIT_S, '_cluster_',cluster_id, '.jpeg')

ggsave(fname,pvol, width=4.5,height=7, dpi=300)


Padj_T_paths=0.05
padj_paths<-Padj_T_paths
pvalueCutoff=1
if (!process_mirnas){
  source('ppmi/RNAseq enrichment.R')
  
}else{
  source('ppmi/miRNA_seq_enrichment.R')
  
}









