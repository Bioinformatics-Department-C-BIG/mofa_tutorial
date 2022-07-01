#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
outdir3<-'Comparisons/plots/mixomics/unsupervised/'
#BiocManager::install("org.Hs.eg.db")

CLL_data$mRNA
library('org.Hs.eg.db')

out_dir<-'Comparisons/plots/'

#### unsupervised is a bit poor in terms of output
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


list.keepX <- list(mRNA = c(16, 16), drug = c(15,16), meth=c(16,16), mut=c(5,5))

# without mutations
list.keepX <- list(mRNA = c(15, 15), mut = c(15,15), meth=c(15,5))

total_comps=15

names(X)


##### TUNING 
spls.result <- block.spls(X, indY=2, ncomp = total_comps, mode='canonical')


# repeated CV tuning of component count
perf.spls.result <- perf(spls.result, validation = 'Mfold',
                        folds = 10, nrepeat = 5) 

plot(spls.result, criterion = 'Q2.total')

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
list.keepX <- list(mRNA = seq(20, 50, 5), mut = seq(20, 50, 5), meth=seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10) 


tune.spls.liver <- tune.spls(X, Y, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 1, folds = 10, # use 10 folds
                             mode = 'canonical', measure = 'cor') 
plot(tune.spls.liver)         # use






### UNSUPERVIED MODEL 
### try to assopciate all others to the drugs! 
block.pls.result <- block.spls(X, indY=2, keepX=list.keepX,ncomp = total_comps)
MyResult.diablo<-block.pls.result

plotIndiv(MyResult.diablo)

plotVar(MyResult.diablo, var.names = c(TRUE, TRUE, TRUE,TRUE),
        legend=TRUE, pch=c(16,16,16,16))


params_str<-paste0('unsup_',paste(unlist(names(X)), collapse='_'), '_', length(names(X)), '_', total_comps)


v_set=c()

for (comp in 1:total_comps){ 
  top_genes_mixomics<-selectVar(MyResult.diablo, comp = comp,list.keepX=10)$mRNA$name[1:list.keepX$mRNA[1]]
  write.csv(top_genes_mixomics, paste0(outdir3, 'top_genes_', total_comps, '_', comp, '.csv'))
  
  top_meth_mixomics<-selectVar(MyResult.diablo, comp = comp,list.keepX=10)$meth$name[1:list.keepX$meth[1]]
  write.csv(top_meth_mixomics, paste0(outdir3, 'top_meth_', total_comps, '_', comp, '.csv'))
  
  
  # concat all 
  v_set <- c(v_set, top_genes_mixomics)
  print(length(v_set))
}
write.csv(v_set, paste0(outdir2, 'all_mutations_top_', total_comps, '.csv'))


#plot the contributions of each feature to each dimension
for (i in 1:5){

  png(paste0(outdir3,'loadings_', total_comps, i, '.png'))
  plotLoadings(block.pls.result, ncomp = i) 
  selectVar(MyResult.diablo, block = 'mRNA', comp = i)$mRNA$name
  dev.off()

  png(paste0(outdir3,'factor_space.png'))
  plotIndiv(block.pls.result, group = Y, comp = 1:2) # plot the samples
  dev.off()

  png(paste0(outdir3,'variable.png'))
  plotVar(block.pls.result, legend =TRUE, comp = 1:2 ) # plot the variables
  dev.off()
  
  
  
  }
















