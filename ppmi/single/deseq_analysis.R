

#install.packages('R.filesets') ; install.packages(c("factoextra", "FactoMineR"))

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
datalist=loadRDS(deseq_file)
ddsSE=datalist[[1]]
vsd=datalist[[2]]
se_filt=datalist[[3]]
deseq2Results=datalist[[4]]

table(se_filt$COHORT_DEFINITION)

# todo join strings
# TODO: Report the number of samples too! 

des<-paste0(as.character(design(ddsSE))[-1])
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

if (run_ma){
  jpeg(paste0(outdir_s, '/MA_plot_results.jpeg'))
  plotMA(deseq2Results)
  dev.off()
}



#, ylim=c(-1,10), xlim=c(0,5))



#### PCA plots

pca.data <- PCA(t(highly_variable_genes_mofa), scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))


fviz_pca_ind(pca.data)

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
  length(ens)
  ens<-gsub('\\..*', '',ens)
  symbols <- mapIds(org.Hs.eg.db, keys = ens,
                    column = c('SYMBOL'), keytype = 'ENSEMBL')
  symbols <- symbols[!is.na(symbols)]
  symbols_ordered <- symbols[match(ens, names(symbols))]
  na_ind<-is.na(symbols_ordered);
  symbols_ordered[na_ind]=ens[na_ind]
  deseq2ResDF$SYMBOL<-symbols_ordered
  
}
rowData(vsd)$SYMBOL=symbols_ordered
is.na(rowData(vsd)$SYMBOL)
#symbols[dup_ind]

outdir_s



# Set a boolean column for significance





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
deseq2ResDF<-mark_signficant(deseq2ResDF, padj_T, log2fol_T)





##### PLOTS of DE genes 
# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3),
                                          oob=squish) + scale_x_log10() +
  geom_hline(yintercept = 0, colour="tomato1", size=2) + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_manual(name="q-value", values=("Significant"="red"),
                      na.value="grey50") + theme_bw()


dir.create(outdir_s)


# Let's add some more detail
dev.off()
limits<-as.numeric(max(abs(deseq2ResDF$log2FoldChange)))
limits
p<-ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-2, 2), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, 
                               linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() +
  geom_density_2d(colour="black", size=0.5)+ 
  geom_hline(yintercept=c(0.15, -0.15),linetype="dashed", color = "red")+
  ylim(-(limits+0.1),(limits+0.1))
  

p

logplot_f<-paste0(outdir_s, '/logfold_padjusted_plot.jpeg')
ggsave(logplot_f,width=8, height=6, res=300 )
ggsave(logplot_f,width=8, height=6 )


#### NORMALIZED COUNTS
top_gene<-rownames(deseq2ResDF[which.min(deseq2ResDF$padj),])
#top_gene

contrasts<-c('Condition', 'EVENT_ID')
dds$Condition<-dds$COHORT
  
# Extract counts for the gene otop2
otop2Counts <- plotCounts(dds, 
                          gene=top_gene, 
                          intgroup=contrasts, 
                          returnData=TRUE)




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
nind<-300
rand_ind<-sample(seq(1:dim(vsd)[2]),nind )
vsd_filt<-vsd[,rand_ind]
vsd_filt=vsd
rowData(vsd)$SYMBOL
deseq2VST <- as.data.frame( assay(vsd_filt))

#### Filter data for visualization 
table(vsd_filt$COHORT)
### First filter data 

dim(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)


deseq2ResDF$padj 
# Keep only the significantly differentiated genes where the fold-change was at least 3
log2fol_T<-0.25
padj_T<-.005

sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= padj_T & abs(deseq2ResDF$log2FoldChange) > log2fol_T,])
orderedSigGenes<-deseq2ResDF[deseq2ResDF$padj <= padj_T  & abs(deseq2ResDF$log2FoldChange) > log2fol_T, ] %>% 
  arrange(padj)
dim(orderedSigGenes)
#  arrange(desc(abs(log2FoldChange)))
ntop=20
sigGenes <- rownames(orderedSigGenes)
#vsd_filt

## Take the first 20 instead of threshold?? 
length(sigGenes)
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

dim(deseq2VST)
#deseq2VST$Gene

#Convert the VST counts to long format for ggplot2
library(reshape2)
library(pheatmap)

deseq2VST_wide <- deseq2VST

deseq2VST_long <- melt(deseq2VST_wide, id.vars=c("Gene"))



head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
#deseq2VST_p <- melt(deseq2VST_long, id.vars=c("Gene"))


graphics.off()
# Make a heatmap ## using the vst
## TODO add annotation

heatmap <- ggplot(deseq2VST_long, aes(x=variable, y=Gene, fill=value)) + geom_raster() + 
  scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
heatmap
library('pheatmap')

df<-vsd_filt$COHORT
assay(vsd_filt)
vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]
#vsd_filt


dim(vsd_filt_genes)
length(vsd_filt_genes$COHORT)

df<-as.data.frame(colData(vsd_filt_genes)[,c("COHORT","SEX", 'NHY')])
df<-as.data.frame(colData(vsd_filt_genes)[,c("COHORT","SEX", 'NHY')])

#colnames(assay(vsd_filt_genes))==vsd_filt_genes$PATNO_EVENT_ID
graphics.off()
fname<-paste0(outdir_s, '/heatmap3', padj_T,'_', log2fol_T ,'.jpeg')
jpeg(fname, width=2000, height=1500, res=200)


my_pheatmap<-pheatmap(assay(vsd_filt_genes), 
        labels_row=as.character(rowData(vsd_filt_genes)$SYMBOL),
         cluster_rows=TRUE, 
         show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df
         )

my_pheatmap

dev.off()

P2<-pheatmap(assay(vsd_filt_genes), 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
P2


fname<-paste0(outdir_s, '/heatmap2.jpeg')


ggsave(paste0(outdir_s, '/heatmap2.jpeg'))




dists <- dist(t(assay(vsd_filt)))
plot(hclust(dists))


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
ggsave(fname)

library('EnhancedVolcano')
pvol+coord_flip()

fname<-paste0(outdir_s, '/EnhancedVolcano_flip.jpeg')
ggsave(fname)











