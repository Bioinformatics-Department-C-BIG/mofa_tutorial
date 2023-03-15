#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
#install.packages('edgeR')
#BiocManager::install('limma')
#### tutorial from: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html


### Problem: our std mean plots are not good : retry with this: 
### https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html

library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
#install.packages('statmod')
library('statmod')
library(gplots)
library(RColorBrewer)
library(sys)
library(sys)



# Run once for proteomics and once for transcriptomics, or change to rn both
# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html

source('bladder_cancer/preprocessing.R')


output_1='bladder_cancer/plots/deseq/'
output_files<-'bladder_cancer/'


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)


prot=FALSE
if (prot){
  seqdata <- read.delim(paste0(dir,'Proteomics_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X2_t_cut) ; rownames(seqdata)<-colnames(X2_t_cut)
  no_name<-which(is.na(seqdata[1]))
  #### Format
  if (length(no_name)>0){
    seqdata<-seqdata[,-no_name]}
  countdata<-seqdata
  
  output_de=paste0(output_1, 'prot')
  ng=ng_p
}else{
  seqdata <- read.delim(paste0(dir,'RNAseq_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X1_t_cut) ; rownames(seqdata)<-colnames(X1_t_cut)
  no_name<-which(is.na(seqdata[1]))
  #### Format
  if (length(no_name)>0){
    seqdata<-seqdata[,-no_name]}
  countdata<-seqdata
  
  
  
  output_de=paste0(output_1, 'genes')
  ng=ng_g
  
}


seqdata
group <- paste(Y_raw$Subtype)
d <- DGEList(counts=seqdata,group=factor(group))

## Calculations will be added to object d !! 
dim(d)
d.full <- d
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum) # total gene counts per sample
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)

### After filtering reset the library sizes 
d$samples$lib.size <- colSums(d$counts)
d$samples
d <- calcNormFactors(d)
d

plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)


#### NORMALIZATION 
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)

png(paste0(output_de, '_ng_', ng ,'_BCV.png'), type='cairo')
plotBCV(d1)
dev.off()



### different tutorial # more complex models
# 2.12 What to do if you have no replicates
#Recall that d1 is the naive method where we only fit a common dispersion.

et12 <- exactTest(d1, pair=c(1,2))
topTags(et12, n=10)
# The total number of differentially expressed genes at FDR< 0:05 is:
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)
## differentially expressed tags from the naive method in d1

de1tags12 <- rownames(d1)[as.logical(de1)] 
plotSmear(et12, de.tags=de1tags12)
abline(h = c(-2, 2), col = "blue")



### glm 

fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1))

write.csv(data.frame(topTags(lrt12, n=20)), paste0(output_de, 'top_genes.csv'))

# https://bioconductor.org/packages/release/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
## Removing heteroscedascity from count data
design <- model.matrix(~0+group)
par(mfrow=c(1,2))
# bfore normalization
v <- voom(d1, design, plot=TRUE)
contr.matrix <- makeContrasts(
 NPS1vsNPS3 = groupNPS1-groupNPS3, 
  levels = colnames(design))
contr.matrix

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

