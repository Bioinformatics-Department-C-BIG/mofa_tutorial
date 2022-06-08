

BiocManager::install('vsn')
BiocManager::install('DEP')

library('vsn')
library('DEP')
## The short answer is that you almost certainly should just log2 
# transform the mass spec values, after adding a constant to avoid 
# zeros, and apply limma and eBayes with trend=TRUE and robust=TRUE. 
# There is no possible advantage in using voom or vsn.


# proteomic analysis vsn 

# Load example
countdata <- seqdata[,-1]

#data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]
#data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Make SummarizedExperiment
columns <- grep("LFQ.", colnames(data_unique))
exp_design <- UbiLength_ExpDesign
exp_design = data.frame(label=Y_raw$Sample,condition=Y_raw$Subtype, replicate=rep(1, 16))
columns=seq(1:16)
se <- make_se(data, columns, exp_design)

# Filter and normalize
filt <- filter_missval(countdata, thr = 0)
norm <- normalize_vsn()

normalized_data<-justvsn(countdata)
meanSdPlot(normalized_data)
fit=vsn2(countdata)
ynorm=predict(fit, countdata)


meanSdPlot(as.matrix(ynorm))

