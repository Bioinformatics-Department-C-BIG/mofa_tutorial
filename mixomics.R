#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("mixOmics")


#browseVignettes("mixOmics")
library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene


MyResult.pca <- pca(X)     # 1 Run the method
plotIndiv(MyResult.pca)    # 2 Plot the samples

plotVar(MyResult.pca)      # 3 Plot the variables

plotIndiv(MyResult.pca, group = liver.toxicity$treatment$Dose.Group, 
          legend = TRUE)



##### N integraiton 

library(mixOmics)
data(breast.TCGA)
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
summary(Y)


list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))


MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)



plotVar(MyResult.diablo, var.names = c(FALSE, FALSE, TRUE),
        legend=TRUE, pch=c(16,16,1))



##### Customize sample plots 

plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2,3),
          title = 'BRCA with DIABLO')

# Only show the names of the proteins 
plotVar(MyResult.diablo, var.names = c(FALSE, FALSE, TRUE),
        legend=TRUE, pch=c(16,16,1))

##### GLoval OVERVIEW of the correlation structure at the componet level 
plotDiablo(MyResult.diablo, ncomp = 1)


### Visualize correlations between variables 
circosPlot(MyResult.diablo, cutoff=0.7)


##### cimDiablo
# Represent the multi-omics molecular signature expression for each sample 



#####

# minimal example with margins improved:
# cimDiablo(MyResult.diablo, margin=c(8,20))

# extended example:
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'), comp = 1, margin=c(8,20), legend.position = "right")



###### Plot loadings: visualizes loading weights of each selected variables on each component
# And each dataset 
#plotLoadings(MyResult.diablo, contrib = "max")
plotLoadings(MyResult.diablo, comp = 2, contrib = "max")


network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')

