#if (!requireNamespace("BiocManager"))
#  install.packages("BiocManager")
#BiocManager::install(c("limma", "edgeR", "Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))


#### tutorial from: https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html


library(edgeR)
library(limma)
library(Glimma)
library(org.Mm.eg.db)
library(gplots)
library(RColorBrewer)
library(NMF)

dir='/Users/efiathieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'

dir='E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'
output_1='bladder_cancer/plots/deseq/'

X1_raw<-read.csv(file = paste0(dir,'RNAseq_BladderCancer.csv' ))
X2_raw<-read.csv(file = paste0(dir,'Proteomics_BladderCancer.csv' ))
Y_raw<-read.csv(file = paste0(dir,'pheno_BladderCancer.csv' ), nrows = 16)

Y_raw$Subtype<-as.factor(Y_raw$Subtype)
Y_raw$Grade<-as.factor(Y_raw$Grade)
Y_raw$TURB.stage<-as.factor(Y_raw$TURB.stage)

prot=FALSE
if (prot){
  seqdata <- read.delim(paste0(dir,'Proteomics_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X1_t_cut)
  countdata<-seqdata
  
  output_de=paste0(output_1, 'prot')
}else{
  seqdata <- read.delim(paste0(dir,'RNAseq_BladderCancer.csv' ), sep=',', stringsAsFactors = FALSE)
  countdata <- seqdata[,-1]
  
  
  seqdata<- t(X2_t_cut)
  countdata<-seqdata
  output_de=paste0(output_1, 'genes')
  
  }


dim(seqdata)


no_name<-which(is.na(seqdata[1]))
#### Format
if (length(no_name)>0){
  seqdata<-seqdata[-no_name,]}
# Remove first  from seqdata


# Look at the output
head(countdata)

rownames(countdata) <- seqdata[,1]


######## convert to dgelist
y <- DGEList(countdata)
y$samples


group <- paste(Y_raw$Subtype,Y_raw$Grade,sep=".")
group <- factor(group)
y$samples$group <- group


group <- paste(Y_raw$Subtype)
group <- factor(group)
y$samples$group <- group

##### 1. Filter lowly expresed genes
# Note that by converting to CPMs we are normalising for 
# the different sequencing depths for each sample.

myCPM <- cpm(countdata)

head(myCPM)
thresh <- myCPM > 0.5
head(thresh)
table(rowSums(thresh))
keep <- rowSums(thresh) >= 2
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

y$samples$lib.size
png(paste0(output_de, '_lib_size.png'))
barplot(y$samples$lib.size,names=colnames(y),las=2)
dev.off()
title("Barplot of library sizes")


######## boxplots
# count data is not normally distributed
# Get log2 counts per million so we can examine count data
logcounts <- cpm(y,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

barplot(y$samples$lib.size,names=colnames(y),las=2)


### multi dimneisional scaleing 
plotMDS(y)
labels <- paste(Y_raw$Sample, Y_raw$Subtype, Y_raw$Grade)


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
glMDSPlot(y, labels=labels, groups=group, folder="mds")


###hierarchical clustering

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(logcounts, 1, var)
head(var_genes)

# Get the gene names for the top 5000 most variable genes
most_var_n=5000

ng_g

select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
select_var_5000 <- names(sort(var_genes, decreasing=TRUE))[1:(length(var_genes)/ng_g)]


# Subset logcounts matrix

highly_variable_lcpm_most_var <- logcounts[select_var_5000,]

if (prot){
  highly_variable_proteins<-highly_variable_lcpm_most_var
  write.csv2(highly_variable_proteins,'highly_variable_proteins_normalized.csv')
}else{
  highly_variable_genes<-highly_variable_lcpm_most_var
  write.csv2(highly_variable_genes,'highly_variable_genes_normalized.csv')
  
}
head(select_var)

dim(as.data.frame(highly_variable_lcpm_most_var))
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)

## Get some nicer colours
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
col.cell <- c("purple","orange")[Y_raw$Subtype]

dev.off()
# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", 
          main="Top 500 most variable genes across samples",
          ColSideColors=col.cell,scale="row", 
          labCol=group)




## normalization
# Apply normalisation to DGEList object
y <- calcNormFactors(y)# normalize
logcounts_norm <- cpm(y,prior.count=2, log=TRUE)  # take cpm

highly_variable_lcpm_norm<-logcounts_norm[select_var_5000,] # and filter!
if (prot){
  highly_variable_proteins_voom<-highly_variable_lcpm_norm
  write.csv2(highly_variable_proteins_voom,'highly_variable_proteins_normalized.csv')
}else{
  highly_variable_genes_voom<-highly_variable_lcpm_norm
  write.csv2(highly_variable_genes_voom,'highly_variable_genes_normalized.csv')
  
}


y$samples

par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 11)
abline(h=0,col="grey")

####### voom transform
par(mfrow=c(1,1))
design <- model.matrix(~ 0 + group)
v <- voom(y,design,plot = TRUE)
par(mfrow=c(1,1))
v <- voom(y,design,plot = TRUE)


par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,main="Unnormalised logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,main="Voom transformed logCPM")
## Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(v$E),col="blue")


#### SAVE
v_most_var<-v$E[select_var_5000,] # and filter!
if (prot){
  highly_variable_proteins<-v_most_var
  write.csv2(highly_variable_proteins_voom,'highly_variable_proteins_normalized_VOOM.csv')
}else{
  highly_variable_genes<-v_most_var
  write.csv2(highly_variable_genes_voom,'highly_variable_genes_normalized_VOOM.csv')
  
}


fit <- lmFit(v)
names(fit)
fit$design
cont.matrix <- makeContrasts(B.NPS3vsNPS1=groupNPS1  - groupNPS3 ,levels=design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)
summary(summa.fit)
library(ggplot2)


### write out results 
# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit.cont,coef=1,status=summa.fit[,"B.NPS3vsNPS1"], 
       values = c(-1, 1), hl.col=c("blue","red"))
dev.off()

# For the volcano plot we have to specify how many of the top genes to highlight.
# We can also specify that we want to plot the gene symbol for the highlighted genes.
# let's highlight the top 100 most DE genes
volcanoplot(fit.cont,coef=1,highlight=20,names=rownames(fit.cont$coefficients),
            main="B.NPS3vsNPS1")
ggsave(paste0(output,'volcano.png'))

