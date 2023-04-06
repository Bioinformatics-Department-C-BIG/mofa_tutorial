

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))

script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir, '/setup_os.R'))
print(script_dir)
library('R.filesets')
library(DESeq2)
library("SummarizedExperiment")
library(data.table)
library(dplyr)

### TODO: Add volcano plot for each time point -DONE
### TODO: add heatmap for all tps tpogether -DONE
#source('ppmi/de')

#load
library("factoextra")
library("FactoMineR")
library('pheatmap')
library('ggplot2')

#### Run DE 


#ddsSE <- DESeqDataSet(se_filt, 
#                      design = ~PATNO)

# TODO: assign the groups 
#dds <- DESeqDataSetFromMatrix(
# countData = assay(se_filt),
#  colData = colData(se_filt),
#  design = ~COHORT, tidy = F
#)

#source(paste0(script_dir, '/config.R'))

### LOAD runs
VISIT='V08'



source(paste0(script_dir, '/config.R'))
source(paste0(script_dir, '/utils.R'))


print(deseq_file)

datalist=loadRDS(deseq_file)
ddsSE=datalist[[1]]
vsd=datalist[[2]]
se_filt=datalist[[3]]
deseq2Results=datalist[[4]]

table(se_filt$COHORT_DEFINITION)

# todo join strings
# TODO: Report the number of samples too! 

des<-gsub(' ', '', paste0(as.character(design(ddsSE))[-1]))
if  (process_mirnas){

  outdir_s<-paste0(outdir_orig, '/single/', param_str_m, des)
  
}else{
  outdir_s<-paste0(outdir_orig, '/single/', param_str_g, des)
  
}
outdir_s
# RUN DIFFERENTIAL EXPRESSION ANALYSIS 


### TODO: save the file so we don't have to fit the model each time!! 
dds<-ddsSE



#rm(deseq2Data)



dir.create(outdir_s)
#deseq2Data<-loadRDS(paste0(outdir_s, '/deseq_results.RDS'))
#### First obtain the single omics significant RNAs 



write.csv(deseq2Results, paste0(outdir_s, '/results.csv'))

deseq2ResDF <- as.data.frame(deseq2Results)


### Up to here output can be used for other reasons
##



# 
run_ma<-FALSE
if (run_ma){
  jpeg(paste0(outdir_s, '/MA_plot_results.jpeg'))
  plotMA(deseq2Results)
  dev.off()
}



#, ylim=c(-1,10), xlim=c(0,5))



#### PCA plots

#pca.data <- PCA(t(highly_variable_genes_mofa), scale.unit = TRUE, graph = FALSE)

#fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))


#fviz_pca_ind(pca.data)

#### Make plots 
# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))
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





mark_signficant<-function(deseq2ResDF, padj_T, log2fol_T){
  ## mark a significant column and write to file
  
  signif_file<-paste0('/significant', padj_T, '_',log2fol_T, '.csv')
  
  deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < padj_T , "Significant", NA)
  deseq2ResDF$sign_lfc <- ifelse(deseq2ResDF$padj <padj_T & abs(deseq2ResDF$log2FoldChange) >log2fol_T , "Significant", NA)
  # Examine this data frame
  # Order the significant to save as a new output file 
  head(deseq2ResDF)
  # LARGER ONE not saved 
  sign_only<-deseq2ResDF[which(deseq2ResDF$sign_lfc=='Significant'),]
  sign_only_ordered<-sign_only[order(sign_only$'padj', decreasing = FALSE),]
  write.csv(sign_only_ordered,paste0(outdir_s,signif_file), row.names = TRUE)
  ### create also a more strict file? 
  return(deseq2ResDF)
}
log2fol_T<-0.25
padj_T<-.005

deseq2ResDF_strict<-mark_signficant(deseq2ResDF, padj_T, log2fol_T)

log2fol_T<-0.1
padj_T<-.05
deseq2ResDF$log2pval<-deseq2ResDF$log2FoldChange*-log10(deseq2ResDF$padj)
deseq2ResDF$abslog2pval<-abs(deseq2ResDF$log2pval)

deseq2ResDF<-mark_signficant(deseq2ResDF, padj_T, log2fol_T)



#hist(-log10(deseq2ResDF$padj))
#hist(deseq2ResDF$log2FoldChange)
#hist(deseq2ResDF$log2pval[abs(deseq2ResDF$log2pval)>0.15])
#
#
#
#min(deseq2ResDF$log2pval, na.rm = TRUE)


###### PLOTS of DE genes 

## Plot the results similar to DEseq2

# Let's add some more detail
dev.off()
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
run_heatmap=TRUE

if (run_heatmap){
  vsd_filt=vsd
  deseq2VST <- as.data.frame( assay(vsd_filt))
  
  #### Filter data for visualization 
  table(vsd_filt$COHORT)
  ### First filter data 
  rownames(vsd)
  rownames(deseq2VST)<-rownames(vsd)
  deseq2VST$Gene <- rownames(deseq2VST)
  head(deseq2VST$Gene)
  
  
  #deseq2ResDF$padj 
  # Keep only the significantly differentiated genes where the fold-change was at least 3
  log2fol_T<-0.15
  padj_T<-.005
  
  sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= padj_T & abs(deseq2ResDF$log2FoldChange) > log2fol_T,])
  deseq2ResDF$Gene<-rownames(deseq2ResDF)
  order_by_metric<-'log2pval'
  order_by_metric<-'abslog2pval'
  
  oSigGenes<-deseq2ResDF[deseq2ResDF$padj <= padj_T  & abs(deseq2ResDF$log2FoldChange) > log2fol_T, ] 
  orderedSigGenes<-oSigGenes[order(-oSigGenes[,order_by_metric]),]

  dim(orderedSigGenes)
  n_sig<-50
  n_sig=dim(orderedSigGenes)[1]
  sigGenes <- orderedSigGenes$Gene[1:n_sig]
  length(sigGenes);head(sigGenes)
  deseq2VST[deseq2VST$Gene %in% sigGenes,]
  deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]
  dim(deseq2VST)
  #deseq2VST$Gene
  
  #Convert the VST counts to long format for ggplot2
  library(reshape2)
  library(pheatmap)
  
  deseq2VST_wide <- deseq2VST
  deseq2VST_long <- melt(deseq2VST_wide, id.vars=c("Gene"))
  
  graphics.off()
  # Make a heatmap ## using the vst
  ## TODO add annotation
  
  
  library('pheatmap')
  
  df<-vsd_filt$COHORT
  assay(vsd_filt)
  vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
  #vsd_filt
  
  
  dim(vsd_filt_genes)
  length(vsd_filt_genes$COHORT)
  
  df<-as.data.frame(colData(vsd_filt_genes)[,c("COHORT","SEX", 'NHY')])
  
  #colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
  graphics.off()
  fname<-paste0(outdir_s, '/heatmap3', '_',padj_T,'_', log2fol_T ,order_by_metric, '_', n_sig,'.jpeg')
  jpeg(fname, width=2000, height=1500, res=200)
  
  if(process_mirnas){
    lab=rownames(rowData(vsd_filt_genes)) }else{
      lab=as.character(rowData(vsd_filt_genes)$SYMBOL)}
  
  
      my_pheatmap<-pheatmap(assay(vsd_filt_genes), 
                            labels_row=lab,
                            cluster_rows=TRUE, 
                            show_rownames=TRUE,
                            cluster_cols=FALSE,
                            annotation_col=df
      )
      
      my_pheatmap
  
  dev.off()
  my_pheatmap
  

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



### Add Volcano plots 
### Compare to their results 




#BiocManager::install('EnhancedVolcano')

library('EnhancedVolcano')
if(process_mirnas){lab=rownames(deseq2ResDF) }else{lab=deseq2ResDF$SYMBOL}


ns<-table(se_filt$COHORT_DEFINITION)
ns<-paste(rownames(ns)[1], ns[1],', ',names(ns)[2],ns[2])
pvol<-EnhancedVolcano(deseq2ResDF,
                lab = lab,
                pCutoff = 10e-6,
                FCcutoff = 0.1,
                x = 'log2FoldChange',
                y = 'pvalue', 
                  
                subtitle = ns   )
pvol
fname<-paste0(outdir_s, '/EnhancedVolcano.jpeg')
ggsave(fname, width=7,height=8)
#
#library('EnhancedVolcano')
#pvol+coord_flip()
#
#fname<-paste0(outdir_s, '/EnhancedVolcano_flip.jpeg')
#ggsave(fname, width=8, height=7)










