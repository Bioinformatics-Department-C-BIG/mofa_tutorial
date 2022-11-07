#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))
#install.packages('edgeR')
#BiocManager::install('limma')
#### tutorial from: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html


library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)

library(gplots)
library(RColorBrewer)
library(sys)
library(sys)


source('bladder_cancer/preprocessing.R')


output_1='bladder_cancer/plots/deseq/'
output_files<-'bladder_cancer/'


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)

prot=FALSE


if (prot){
  lowest_thresh<-10
  params.remove_low_threshold<-15
  params.ng<-2
}else{
  lowest_thresh<-10
  params.remove_low_threshold<-10
  params.ng<-5
  

}



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
  
  }


dim(seqdata)



# Remove first  from seqdata


# Look at the output
head(countdata)



######## convert to dgelist
y <- DGEList(countdata)
y$samples


group_all <- paste(Y_raw$Subtype,Y_raw$Grade,sep=".")
group <- paste(Y_raw$Subtype)

y$samples$group <- group

group_all <- factor(group)


group <- factor(group)
y$samples$group <- group



y$samples$Age <- paste(Y_raw$Age)
y$samples$Sex <- factor(paste(Y_raw$Sex))

##### 1. Filter lowly expresed genes
# Note that by converting to CPMs we are normalising for 
# the different sequencing depths for each sample.

myCPM <- cpm(countdata)

head(myCPM)



thresh <- myCPM > lowest_thresh
head(thresh)
table(rowSums(thresh))

keep <- rowSums(thresh) >= params.remove_low_threshold
summary(keep)


# Let's have a look and see whether our threshold of 0.5 does indeed correspond to a count of about 10-15
# We will look at the first sample
# Let us limit the x and y-axis so we can actually look to see what is happening at the smaller counts
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,30))
# Add a vertical line at 0.5 CPM
abline(v=0.5)


# Filter 
y <- y[keep, keep.lib.sizes=FALSE]

##### QC: Library size

png(paste0(output_de, '_lib_size.png'), type='cairo')
barplot(y$samples$lib.size,names=colnames(y),las=2)
dev.off()
title("Barplot of library sizes")


######## boxplots
# count data is not normally distributed
# Get log2 counts per million so we can examine count data

logcounts <- cpm(y,log=TRUE)
if (prot){
  logcounts<-log(y$counts)
}
hist(logcounts)

# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

barplot(y$samples$lib.size,names=colnames(y),las=2)


### multi dimneisional scaleing 
plotMDS(y)
labels <- paste(Y_raw$Sample, Y_raw$Subtype, Y_raw$Grade)

NROW(y)
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,2))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(Y_raw$Subtype)

## Let's choose purple for basal and orange for luminal
col.cell <- c("purple","orange")[as.factor(Y_raw$Subtype)]
data.frame(Y_raw$Subtype,col.cell)
# Redo the MDS with cell type colouring
plotMDS(y,col=col.cell)
# Let's add a legend to the plot so we know which colours correspond to which cell type
legend("topleft",fill=c("purple","orange"),legend=levels(Y_raw$Subtype))
# Add a title
title("Cell type")

# Similarly for status
levels(Y_raw$TURB.stage)
col.status <- c("blue","red","black")[Y_raw$TURB.stage]
plotMDS(y,col=col.status)
legend("topleft",fill=c("blue","red","black"),legend=levels(Y_raw$TURB.stage),cex=0.8)
title("Status")
col.status

labels <- paste(Y_raw$TURB.stage, Y_raw$Subtype,  Y_raw$Grade)
#glMDSPlot(y, labels=labels, groups=group, folder="mds")


###hierarchical clustering


# We estimate the variance for each row in the logcounts matrix

all_variances <- apply(logcounts, 1, var)

select_most_variable <- names(sort(all_variances, decreasing=TRUE))[1:(length(all_variances)/  params.ng)]

highly_variable_mofa<-logcounts[select_most_variable,]


#highly_variable_mofa<-highly_variable_mofa[highly_variable_mofa>2,]
hist(highly_variable_mofa)
  

hist(highly_variable_mofa)

if (prot){
  highly_variable_proteins_mofa<-highly_variable_mofa
  write.csv(highly_variable_proteins_mofa,'highly_variable_proteins_mofa.csv')
}else{
  highly_variable_genes_mofa<-highly_variable_mofa
  write.csv(highly_variable_genes_mofa,'highly_variable_genes_mofa.csv')
}

hist(highly_variable_proteins_mofa)
hist(highly_variable_genes_mofa)
NROW(highly_variable_genes_mofa)

