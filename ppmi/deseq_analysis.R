

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))
source(paste0(script_dir,'/deseq_analysis_setup.R'))

process_mirnas
write.csv(deseq2Results, paste0(outdir_s, '/results.csv'))

deseq2ResDF <- as.data.frame(deseq2Results)
deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
deseq2ResDF$abslog2pval<-abs(deseq2ResDF$log2pval)

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


#symbols[dup_ind]

outdir_s

deseq2ResDF$SYMBOL




log2fol_T<-0.25
padj_T<-.005

deseq2ResDF_strict<-mark_signficant(deseq2ResDF, padj_T, log2fol_T)


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
  for (most_var in c(0.05, 0.1,0.15,0.2,0.25,0.3,  0.9,0.75,0.5)){
    

  param_str_tmp<-paste0(prefix, VISIT_S, '_',most_var ,'_', min.count, '_coh_', sel_coh_s, '_'  )
  highly_variable_outfile<-paste0(output_files, param_str_tmp,'_highly_variable_genes_mofa.csv')
  highly_variable_sign_outfile<-paste0(output_files, param_str_tmp,'_highly_variable_genes_mofa_signif.csv')

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
deseq2ResDF<-mark_signficant(deseq2ResDF, padj_T_overall, log2fol_T_overall, outdir_single = outdir_s)

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


# Plot the data using ggplot2
# this is more useful for paired datasets !! 
# TODO: otherwise do a box plot 
#p<-ggplot(otop2Counts, aes(x=Condition, y=count, colour=Sample, 
#                           group=Sample)) + geom_point() + geom_line() +
#  theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) +
#  guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")
#
#
#p
#
#p<-ggplot(otop2Counts, aes(x=Condition, y=count, colour=Sample, 
#                           group=Sample)) + geom_point() + geom_line() +
#  theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) +
#  guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")
#
#
#p

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
  oSigGenes<-oSigGenes[oSigGenes$baseMean>mean_expr_T,];dim(oSigGenes)
  
  ### Order using the order_by_metric 
  orderedSigGenes<-oSigGenes[order(-oSigGenes[,order_by_metric]),]
  
  n_sig_f='all'
  n_sig_f=30
  
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
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]

  
  ### Add the annotations 
  df<-as.data.frame(colData(vsd_filt_genes)[,c("COHORT", "SEX", 'AGE', 'NHY')])
  
  ## HEATMAP OPTIONS 
  cluster_cols=TRUE
    #colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
  fname<-paste0(outdir_s, '/heatmap3', '_',padj_T_hm,'_', log2fol_T_hm ,order_by_metric, 'high_var_' ,
                filter_highly_var,    '_', most_var_t, '_',  n_sig_f, cluster_cols, '.jpeg')
  
  #ARRANGE
  df_ord<-df[order(df$COHORT),]
  hm<-assay(vsd_filt_genes)
  hm_ord<-hm[,order(df$COHORT)]
  
  ### SCALE!! 
  
  hm_scaled <- as.matrix(hm_ord) %>% t() %>% scale() %>% t()
  dim(hm_ord)
  #jpeg(fname, width=2000, height=1500, res=200)
  graphics.off()
  library(ggplot2)
  if(process_mirnas){
    lab=rownames(rowData(vsd_filt_genes)) }else{
      lab=as.character(rowData(vsd_filt_genes)$SYMBOL)}
  
      #jpeg(fname, width=10*100, height=7*100, res=300)
      my_pheatmap<-pheatmap(hm_scaled, 
                            labels_row=lab,
                            cluster_rows=TRUE, 
                            show_rownames=TRUE,
                            #scale='row', 
                            cluster_cols=cluster_cols,
                            annotation_col=df_ord, 
                              clustering_method = 'ward.D2'
      )
      
      
     show(my_pheatmap)
      my_pheatmap
     # ggsave(fname, width=10, height=7)
      
      ggsave(fname,plot=my_pheatmap, width=7, height=7, dpi=300)


  #P2<-pheatmap(assay(vsd_filt_genes), 
  #         cluster_rows=FALSE, 
  #         show_rownames=TRUE,
  #         cluster_cols=TRUE, annotation_col=df)
  #P2
  
  
  #fname<-paste0(outdir_s, '/heatmap2.jpeg')
  #ggsave(fname, width=8, height=8)
  
  #dev.off()
  
  
  
  #dists <- dist(t(assay(vsd_filt)))
  #plot(hclust(dists))
}

my_pheatmap

### Add Volcano plots 
### Compare to their results 



deseq2ResDF$SYMBOL

library('EnhancedVolcano')
if(process_mirnas){lab=rownames(deseq2ResDF) }else{lab=deseq2ResDF$SYMBOL}

mfc<-max(abs(deseq2ResDF$log2FoldChange))
pmax<-max(-log10(deseq2ResDF$padj), na.rm = TRUE)
pmax
xlim = c(-mfc-0.2,mfc+0.2)
#ylim = c(0,pmax+1)

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

#library(gridExtra)
#library(grid)
##grid.arrange(pvol, p2,
#grid.arrange(pvol,
#             ncol=1,
#             top = textGrob('EnhancedVolcano',
#                            just = c('center'),
#                            gp = gpar(fontsize = 32))
#             )
#
#
#
#fname<-paste0(outdir_s, '/EnhancedVolcano.jpeg')
#ggsave(fname, width=9,height=8)
##
#library('EnhancedVolcano')
#pvol+coord_flip()
#
#fname<-paste0(outdir_s, '/EnhancedVolcano_flip.jpeg')
#ggsave(fname, width=8, height=7)
Padj_T_paths=0.05
padj_paths<-Padj_T_paths
pvalueCutoff=1
if (!process_mirnas){
  source('ppmi/RNAseq enrichment.R')
  
}









