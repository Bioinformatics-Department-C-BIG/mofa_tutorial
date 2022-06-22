#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("mixOmics")
CLL_data$mRNA
library('org.Hs.eg.db')

out_dir<-'Comparisons/plots/'

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
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
summary(Y)

X <- list(mRNA = t(CLL_data$mRNA), 
            drug = t(CLL_data$Drugs))
Y <- CLL_metadata$IGHV

# FILTER MISSING VALUES 
tokeep<-which(!is.na(CLL_metadata$IGHV))



X <- list(mRNA = t(CLL_data$mRNA[,tokeep]), 
          drug = t(CLL_data$Drugs[,tokeep]), 
          meth=t(CLL_data$Methylation[,tokeep])
          #mut=t(CLL_data$Mutations[,tokeep])
          )

ensemblsIDS<-colnames(X$mRNA)
symbols <- mapIds(org.Hs.eg.db, keys = ensemblsIDS, keytype = "ENSEMBL", column="SYMBOL")


not_na<-which(!is.na(symbols))
colnames(X$mRNA)[not_na]<-symbols[not_na]
Y<-as.factor(CLL_metadata$IGHV[tokeep])

NROW(unique(X$mRNA))
# drop a duplicate
ind_to_drop<-which(duplicated(colnames(X$mRNA)))
new<-X$mRNA[,-c(ind_to_drop)]

duplicated(colnames(new))
X$mRNA<-new
#list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))

MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX)
plotIndiv(MyResult.diablo)



plotVar(MyResult.diablo, var.names = c(TRUE, TRUE, TRUE),
        legend=TRUE, pch=c(16,16,16))


# print the first component 
selectVar(MyResult.diablo, block = 'mRNA', comp = 1)$mRNA$name
selectVar(MyResult.diablo, block = 'meth', comp = 1)$meth$name


# todo WHAT can we do with the methylation results now? 

##### Customize sample plots 

plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2),
          title = 'CLL with DIABLO')



##### GLoval OVERVIEW of the correlation structure at the componet level 
plotDiablo(MyResult.diablo, ncomp = 1)

plotArrow(MyResult.diablo, ind.names = FALSE, legend = TRUE, title = 'DIABLO')


### Visualize correlations between variables 
p<-circosPlot(MyResult.diablo, cutoff=0.9)
save(file= 'circos.png' )


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
comp=1
p<-plotLoadings(MyResult.diablo, comp = comp, contrib = "max")
ggsave(paste0(out_dir,'top_weights', comp, '.png'), width = 4, height = 4, dpi=500 )
dev.off()



P<-plotLoadings(MyResult.diablo, comp = 2, contrib = "max")
ggsave('top_weights.png', width = 4, height = 4, dpi=500 )




?network

network(MyResult.diablo, blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')

