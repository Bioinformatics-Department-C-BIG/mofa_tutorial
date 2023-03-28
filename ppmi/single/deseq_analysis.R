#install.packages('R.filesets')
library('R.filesets')


### TODO: Add volcano plot for each time point
### TODO: add heatmap for all tps tpogether 
#source('ppmi/de')


#intall
#install.packages(c("factoextra", "FactoMineR"))

#load
library("factoextra")
library("FactoMineR")
library('pheatmap')

#### Run DE 


#ddsSE <- DESeqDataSet(se_filt, 
#                      design = ~PATNO)



# TODO: assign the groups 
#dds <- DESeqDataSetFromMatrix(
# countData = assay(se_filt),
#  colData = colData(se_filt),
#  design = ~COHORT, tidy = F
#)

highly_variable_genes_mofa
se_filt
ddsSE
dds
# todo join strings
deseq2Data <- DESeq(ddsSE)

des<-paste0(as.character(design(ddsSE))[-1])
sel_coh
if  (process_mirnas){
  outdir_s<-paste0(outdir_orig, '/single/', param_str_m, 'visits_', VISIT_S, '_coh_', sel_coh_s, '_',des)
  
}else{
  outdir_s<-paste0(outdir_orig, '/single/', param_str_g, 'visits_', VISIT_S, '_coh_', sel_coh_s, '_',des)
  
}
outdir_s
# RUN DIFFERENTIAL EXPRESSION ANALYSIS 


### TODO: save the file so we don't have to fit the model each time!! 
dds<-ddsSE



#rm(deseq2Data)



dir.create(outdir_s)
#deseq2Data<-loadRDS(paste0(outdir_s, '/deseq_results.RDS'))
#### First obtain the single omics significant RNAs 

print(deseq2Data$EVENT_ID)


deseq2Results <- results(deseq2Data)
write.csv(deseq2Results, paste0(outdir_s, '/results.csv'))


summary(deseq2Results)

# 
dev.off()

jpeg(paste0(outdir_s, '/MA_plot_results.jpeg'))
plotMA(deseq2Results)
dev.off()
#, ylim=c(-1,10), xlim=c(0,5))



#### PCA plots

pca.data <- PCA(t(highly_variable_genes_mofa), scale.unit = TRUE, graph = FALSE)

fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))

dev.off()
fviz_pca_ind(pca.data)

dev.off()
#### Make plots 
# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))
library(ggplot2)

library(scales) # needed for oob parameter
library(viridis)
#install.packages('viridis')
#install.packages('gridExtra')


# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)

# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .05 , "Significant", NA)

deseq2ResDF$significant 

# Examine this data frame
# Order the significant to save as a new output file 
head(deseq2ResDF)
sign_only<-deseq2ResDF[which(deseq2ResDF$significant=='Significant'),]
sign_only_ordered<-sign_only[order(sign_only$'padj', decreasing = FALSE),]

ens_ids<-gsub('\\..*', '', rownames(sign_only_ordered))
rownames(sign_only_ordered)<-ens_ids

write.csv(sign_only_ordered,paste0(outdir_s, '/significant.csv'), row.names = TRUE)




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
max(deseq2ResDF$log2FoldChange)

p<-ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-2, 2), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, 
                               linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() +
  geom_density_2d(colour="black", size=0.5)


p
logplot_f<-paste0(outdir_s, '/logfold_padjusted_plot.jpeg')
ggsave(logplot_f)

#### NORMALIZED COUNTS
top_gene<-rownames(deseq2ResDF[which.min(deseq2ResDF$padj),])
top_gene

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
p<-ggplot(otop2Counts, aes(x=Condition, y=count, colour=Sample, 
                           group=Sample)) + geom_point() + geom_line() +
  theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) +
  guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")


p

p<-ggplot(otop2Counts, aes(x=Condition, y=count, colour=Sample, 
                           group=Sample)) + geom_point() + geom_line() +
  theme_bw() + theme(axis.text.x=element_text(angle=15, hjust=1)) +
  guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")


p







##### 
# Compute normalization factors and vst 



dds <- estimateSizeFactors(dds,)
# Variance stabilization transformation
# This uses the size factors estimated before 
# TODO: you can run VST using a saved dispersion function
vsd <- varianceStabilizingTransformation(dds)

rand_ind<-sample(seq(1:dim(vsd)[2]), 100)
vsd_filt<-vsd[,rand_ind]

deseq2VST <- as.data.frame( assay(vsd_filt))

#### Filter data for visualization 
table(vsd_filt$COHORT)
### First filter data 

dim(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)


deseq2ResDF$padj 
# Keep only the significantly differentiated genes where the fold-change was at least 3
log2fol_T<-0.35
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > log2fol_T,])
#sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05, ])

length(sigGenes)
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

dim(deseq2VST)
deseq2VST$Gene

#Convert the VST counts to long format for ggplot2
library(reshape2)
library(pheatmap)

deseq2VST_wide <- deseq2VST

deseq2VST_long <- melt(deseq2VST_wide, id.vars=c("Gene"))



head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
#deseq2VST_p <- melt(deseq2VST_long, id.vars=c("Gene"))


dim(deseq2VST_p)
graphics.off()
# Make a heatmap ## using the vst
## TODO add annotation

heatmap <- ggplot(deseq2VST_long, aes(x=variable, y=Gene, fill=value)) + geom_raster() + 
  scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
heatmap
#library('pheatmap')

df<-vsd_filt$COHORT
assay(vsd_filt)
vsd_filt_genes <- vsd_filt[rownames(vsd_filt) %in% sigGenes,]

dim(vsd_filt_genes)
length(vsd_filt_genes$COHORT)

df<-as.data.frame(colData(vsd_filt_genes)[,c("COHORT","SEX")])

pheatmap(assay(vsd_filt_genes), cluster_rows=FALSE, 
         show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df )



ggsave(paste0(outdir_s, '/heatmap2.jpeg'))



dists <- dist(t(assay(vsd_filt)))
plot(hclust(dists))

