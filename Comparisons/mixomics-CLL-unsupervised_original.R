#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db")

CLL_data$mRNA
library('org.Hs.eg.db')

out_dir<-'Comparisons/plots/'


outdir2<-'Comparisons/plots/mixomics/'
#browseVignettes("mixOmics")

library(mixOmics)


data(liver.toxicity)
CLL_data



MyResult.pca <- pca(X)     # 1 Run the method
plotIndiv(MyResult.pca)    # 2 Plot the samples

plotVar(MyResult.pca)      # 3 Plot the variables

plotIndiv(MyResult.pca, group = liver.toxicity$treatment$Dose.Group, 
          legend = TRUE)



##### N integraiton 

library(mixOmics)
data(breast.TCGA)
# DIABLO SUPERVISED
# extract training data and name each data frame

# install.packages("gtools") ## Uncomment if not already installed
library(gtools)
mixedorder( CLL_metadata$sample)

X <- list(mRNA = t(CLL_data$mRNA), 
            drug = t(CLL_data$Drugs))
Y <- CLL_metadata$IGHV[]




# FILTER MISSING VALUES 

rownames(CLL_data$mRNA)

new_order<-order(match(CLL_metadata$sample, colnames(CLL_data$mRNA)))
CLL_metadata_new_order<-CLL_metadata[new_order, ]
tokeep<-which(!is.na(CLL_metadata_new_order$IGHV))

X <- list(mRNA = t(CLL_data$mRNA[,tokeep]), 
          drug = t(CLL_data$Drugs[,tokeep]), 
          meth=t(CLL_data$Methylation[,tokeep]),
          mut=t(CLL_data$Mutations[,tokeep])
          )

ensemblsIDS<-colnames(X$mRNA)
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")


not_na<-which(!is.na(symbols))
colnames(X$mRNA)[not_na]<-symbols[not_na]
Y<-as.factor(CLL_metadata_new_order$IGHV[tokeep])

NROW(unique(X$mRNA))
# drop a duplicate
ind_to_drop<-which(duplicated(colnames(X$mRNA)))
new<-X$mRNA[,-c(ind_to_drop)]
duplicated(colnames(new))
X$mRNA<-new


#list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))

list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), drug=c(15,15))

list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), mut=c(15,15))

list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), mut=c(15,5))
### SANITY CHECKS 
IGHV_stats<-data.frame(x=X$mut[,'IGHV'],y=Y)
#IGHV_stats<-data.frame(x=CLL_metadata$IGHV,
#y=CLL_data$Mutations['IGHV',])

CLL_data$Mutations['IGHV',]

total_comps=5
MyResult.diablo <- block.spls(X, keepX=list.keepX,ncomp = total_comps)
