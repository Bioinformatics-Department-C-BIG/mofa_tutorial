#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
outdir3<-'Comparisons/plots/mixomics/unsupervised/'
#BiocManager::install("org.Hs.eg.db")

CLL_data$mRNA
library('org.Hs.eg.db')

out_dir<-'Comparisons/plots/'


outdir2<-'Comparisons/plots/mixomics/'
#browseVignettes("mixOmics")

library(mixOmics)
# http://mixomics.org/methods/multiblock-spls/

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
#remove missing Ys because it is supervised 
tokeep<-which(!is.na(CLL_metadata_new_order$IGHV))


X <- list(mRNA = t(CLL_data$mRNA[,tokeep]), 
          drug = t(CLL_data$Drugs[,tokeep]), 
          meth=t(CLL_data$Methylation[,tokeep]),
          mut=t(CLL_data$Mutations[,tokeep]))

ensemblsIDS<-colnames(X$mRNA)
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")


not_na<-which(!is.na(symbols))
colnames(X$mRNA)[not_na]<-symbols[not_na]
Y<-as.factor(CLL_metadata_new_order$IGHV[tokeep])

NROW(unique(X$mRNA))
# drop a duplicate
ind_to_drop<-which(duplicated(colnames(X$mRNA)))
new<-X$mRNA[,-c(ind_to_drop)]
# remove duplicate names from X genes
X$mRNA<-new

# Remove NAs in Y :


#list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))

list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), drug=c(15,15))


list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), mut=c(15,15))

# without mutations
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))

total_comps=2

names(X)

### UNSUPERVIED MODEL 
### try to assopciate all others to the drugs! 
block.pls.result <- block.spls(X, indY=2, keepX=list.keepX,ncomp = total_comps)
plotIndiv(MyResult.diablo)

plotVar(MyResult.diablo, var.names = c(TRUE, TRUE, TRUE),
        legend=TRUE, pch=c(16,16,16))


params_str<-paste0('unsup_',paste(unlist(names(X)), collapse='_'), '_', length(names(X)), '_', total_comps)


#plot the contributions of each feature to each dimension
png(paste0(outdir3,'loadings.png'))
plotLoadings(block.pls.result, ncomp = 1) 
dev.off()

png(paste0(outdir3,'factor_space.png'))
plotIndiv(block.pls.result, group = Y) # plot the samples
dev.off()

png(paste0(outdir3,'variable.png'))
plotVar(block.pls.result, legend = TRUE) # plot the variables
dev.off()
