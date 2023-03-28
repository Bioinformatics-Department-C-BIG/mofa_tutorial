
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)
library(DESeq2)
library("vsn")

library("SummarizedExperiment")

library(data.table)
library(dplyr)
## Output directory

output_1='ppmi/output/'
output_files_orig<-'ppmi/output/'

output_de=paste0(output_1, 'gene')
source('bladder_cancer/preprocessing.R')
# TODO: move the pre-processing script to utils


##### Load required data 
# TODO: input all the visits 

MIN_COUNT_G=100
MIN_COUNT_M=10
TOP_GN=0.10
TOP_MN=0.50

VISIT='BL'
VISIT='V04'
VISIT='V08'

g_params<-paste0(VISIT, '_', TOP_GN, '_', MIN_COUNT_G, '_')
m_params<-paste0(VISIT, '_', TOP_MN, '_', MIN_COUNT_M, '_') 




metadata_output<-paste0(output_files, 'combined.csv')
combined<-read.csv2(metadata_output)

#### Remove low expression 
process_mirnas<-FALSE
if (process_mirnas){
  VISIT='BL'
   mirnas_file<-paste0(output_files, 'mirnas_',VISIT,  '.csv')
   mirnas_BL<-as.data.frame( as.matrix(fread(mirnas_file, header=TRUE), rownames=1))
   
   VISIT='V08'
   mirnas_file<-paste0(output_files, 'mirnas_',VISIT,  '.csv')
   mirnas_V08<-as.data.frame(as.matrix(fread(mirnas_file, header=TRUE), rownames=1))
   
  raw_counts<-as.data.frame(mirnas_BL)

  # if we filter too much we get normalization problems 
  min.count=MIN_COUNT_M
  most_var=TOP_MN
  param_str_m<-paste0('mirnas_', m_params)
  highly_variable_outfile<-paste0(output_files, param_str_m,'_highly_variable_genes_mofa.csv')
  
  
}else{
  VISIT='BL'
  rnas_file<-paste0(output_files, 'rnas_', VISIT, '.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
  df_T1<-as.data.frame(rnas_BL)
  
  VISIT='V06'
  rnas_file<-paste0(output_files, 'rnas_', VISIT, '.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
  df_T2<-as.data.frame(rnas_BL)
  
  VISIT='V08'
  rnas_file<-paste0(output_files, 'rnas_', VISIT, '.csv')
  rnas_BL<-as.matrix(fread(rnas_file, header=TRUE), rownames=1)
  df_T3<-as.data.frame(rnas_BL)
  

  
    # this is defined later but filter here if possible to speed up
  # TODO: fix and input common samples as a parameter
# raw_counts<-raw_counts %>% select(common_samples)

  min.count=MIN_COUNT_G
  most_var=TOP_GN
  param_str_g<-paste0('rnas_', g_params )
  highly_variable_outfile<-paste0(output_files, param_str_g,'_highly_variable_genes_mofa.csv')
  outdir_s<-paste0(outdir_orig, '/single/', param_str_g, VISIT)
  
  
  
}



#### TODO: 
# BIND ALL VISITS AND TIME POINTS TOO!! 



# ONLY BIND WHATS common
common_1_2<-intersect(colnames(df_T1), colnames(df_T2))
df_T1<-df_T1 %>% select(all_of(common_1_2))
df_T2<-df_T1 %>% select(all_of(common_1_2))


tps_merged<-merge(df_T1,df_T2, by='row.names', suffix=c('_BL', '_V08'))
rownames(tps_merged)<-tps_merged$Row.names
tps_merged$Row.names<-NULL


Sample<-as.factor(c(colnames(df_T1),colnames(df_T2 )))



time<-c( rep('BL', length(df_T1)), rep('V08', length(df_T2) )) 



#raw_counts<-mutate_all(raw_counts, as.numeric)
raw_counts<-tps_merged

hist(log10(as.matrix(raw_counts)))
dev.off()
## filterbyExpr takes cpm so remove from there 
#cpm_data<-cpm(countdata, log=FALSE)



##### Define

### TODO: Question: Should I input everything into the matrix to normalize? 
### And then filter 

### batch effect and normalization 
# Create a separate matrix with counts only
counts_only <- raw_counts
# Include batch information if there is any
#sample_info$Batch <- as.factor(sample_info$Batch)





### VST? The variance stabilizing and rlog transformations are provided for applications other than differential testing, for example clustering of samples or other machine learning applications. For differential testing we recommend the DESeq function applied to raw counts as outlined above.



#### NUMBER 2 COMPARE WITH CONTROLS 
##### separate cohorts and create filters 

# select patients 
dim(rnas_BL)


raw_counts<-as.data.frame(rnas_BL)

#raw_counts<-as.data.frame(tps_merged)

idx <- edgeR::filterByExpr(raw_counts,min.count=min.count)

length(which(idx))
raw_counts <- as.matrix(raw_counts[idx, ])
dim(raw_counts)

df<-raw_counts


combined_bl<-combined[combined$EVENT_ID==VISIT,]
combined_bl<-combined_bl[combined_bl$COHORT %in% c(1,2),]
dim(df); dim(combined_bl)
common_in_meta<-unique(intersect( colnames(df), combined_bl$PATNO ))
length(common_in_meta)
mirnas_1_metadata<-combined_bl[match(common_in_meta, combined_bl$PATNO),]

df_filt<-df[, match(common_in_meta, colnames(df))] # select only columns with common patients 


dim(mirnas_1_metadata)
dim(df_filt)
combined_bl$PATNO

counts_only<-df_filt
Sample<-as.factor(c(colnames(df_filt)))
Condition<-as.factor(mirnas_1_metadata$COHORT)
Condition


#### common pipeline again
sample_info<-DataFrame(Sample=Sample, Condition=Condition)

# TODO: assign the groups 
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_only),
  colData = sample_info,
  design = ~ Condition, tidy = F
)


#### Run DE 
deseq2Data <- DESeq(dds)



deseq2Results <- results(deseq2Data)

summary(deseq2Results)


plotMA(deseq2Results, ylim=c(-1,10), xlim=c(0,5))








# Load libraries
# install.packages(c("ggplot2", "scales", "viridis"))
library(ggplot2)
library(scales) # needed for oob parameter
library(viridis)

# Coerce to a data frame
deseq2ResDF <- as.data.frame(deseq2Results)



# Set a boolean column for significance
deseq2ResDF$significant <- ifelse(deseq2ResDF$padj < .05 , "Significant", NA)

deseq2ResDF$significant 

# Examine this data frame
head(deseq2ResDF)
sign_only<-deseq2ResDF[which(deseq2ResDF$significant=='Significant'),]
sign_only_ordered<-sign_only[order(sign_only$'padj', decreasing = FALSE),]

write.csv(sign_only_ordered,paste0(outdir_s, '_significant.csv'), row.names = FALSE)


# Plot the results similar to DEseq2
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=significant)) +
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3),
                                        oob=squish) + scale_x_log10() +
  geom_hline(yintercept = 0, colour="tomato1", size=2) + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_manual(name="q-value", values=("Significant"="red"),
                      na.value="grey50") + theme_bw()

# Let's add some more detail
ggplot(deseq2ResDF, aes(baseMean, log2FoldChange, colour=padj)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) +
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, 
                               linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change") +
  scale_colour_viridis(direction=-1, trans='sqrt') + theme_bw() +
  geom_density_2d(colour="black", size=0.05)





#### NORMALIZED COUNTS
top_gene<-rownames(deseq2ResDF[which.min(deseq2ResDF$padj),])
# Extract counts for the gene otop2
otop2Counts <- plotCounts(dds, 
                          gene=top_gene, 
                          intgroup=c('Condition'), 
                          returnData=TRUE)




# Plot the data using ggplot2
# this is more useful for paired datasets !! 
# TODO: otherwise do a box plot 
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



deseq2VST <- as.data.frame( assay(vsd))
deseq2VST$Gene <- rownames(deseq2VST)
head(deseq2VST)


deseq2ResDF$padj 
# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(deseq2ResDF[deseq2ResDF$padj <= .05 & abs(deseq2ResDF$log2FoldChange) > 1,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]



#Convert the VST counts to long format for ggplot2
library(reshape2)

deseq2VST_wide <- deseq2VST
deseq2VST_long <- melt(deseq2VST, id.vars=c("Gene"))


head(deseq2VST_wide)
head(deseq2VST_long)

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap
heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + geom_raster() + 
  scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), 
        axis.ticks.y=element_blank())
heatmap



dists <- dist(t(assay(vsd)))
plot(hclust(dists))



vsd_mat <- assay(vsd)
colnames(vsd_mat)<-vsd$Sample

meanSdPlot(vsd)
dev.off()


DESeq(assay(vsd))




##### Checks
# Check the effect of vst before and after
par(mfrow=c(1,3))
# Check distributions of samples using boxplots
boxplot(log2(assay(dds)[,1:30]), xlab="", ylab="Log2 counts ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(log10(assay(dds))),col="blue")
title("Boxplots of logCPMs (unnormalised)")
boxplot(log10(raw_counts)[,1:30], xlab="", ylab="Log10 counts ",las=2)

# Check distributions of samples using boxplots
boxplot(vsd_mat[,1:30], xlab="", ylab="vst(counts) ",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(vsd_mat),col="blue")
title("Boxplots of logCPMs (after vst)")

## ASSESS BATCH
### Assess batch effect
#plotPCA(vsd, "sample") + labs(color='sample') + ggtitle("Batch effect") 
#dev.off()

##### Store the most variable genes only for MOFA 
# Select most variable genes
highly_variable_genes_mofa<-selectMostVariable(vsd_mat, most_var)
dim(highly_variable_genes_mofa)


boxplot(highly_variable_genes_mofa[,1:30], xlab="", ylab="vst(counts) ",las=2)


write.csv(highly_variable_genes_mofa, highly_variable_outfile, col.names = TRUE)


dim(highly_variable_genes_mofa)

# Check that the distribution is approximately normal
dev.off()

hist(highly_variable_genes_mofa)




### SANITY CHECK: Just plot one gene before and after to ensure the mapping looks correct 
df=highly_variable_genes_mofa
idx=30
plot(df[idx,1:20]) 
gname<-rownames(df)[idx]
title(gname)
df=raw_counts
plot(df[gname, 1:20])
title(gname)

