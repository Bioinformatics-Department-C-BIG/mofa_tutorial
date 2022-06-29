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
#remove missing Ys because it is supervised 
tokeep<-which(!is.na(CLL_metadata_new_order$IGHV))


X <- list(mRNA = t(CLL_data$mRNA[,tokeep]), 
          drug = t(CLL_data$Drugs[,tokeep]), 
          meth=t(CLL_data$Methylation[,tokeep]))

          #mut=t(CLL_data$Mutations[,tokeep]))



# renaming
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


new_drug_names<-rep(drug_comps$Compound, times=1, each=5)
length(new_drug_names)
rep(seq(1:5), times=length(drug_comps) )
length(X$drug[1,])
# Remove NAs in Y :


#list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))

list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), drug=c(15,15))


list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), mut=c(15,15))

# without mutations
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))


names(X)
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX,ncomp = total_comps)
plotIndiv(MyResult.diablo)

basic.diablo.model<- block.splsda(X, Y, ncomp = total_comps)

##### TUNING 
# 1. tune the number of components 


perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 4, nrepeat = 4) 

png(paste0(outdir2, '_tuning_','.png'))
plot(perf.diablo) # plot output of tuning

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 

ncomp=2
total_comps=ncomp
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 



### Tune the number of features 
# set grid of values for each component to test
test.keepX = list (mRNA = c(6:8, seq(10, 18, 5), seq(20,30,6)), 
                   meth = c(6:8, seq(10, 18, 5), seq(20,30,6)),
                   drug = c(6:8, seq(10, 18, 5), seq(20,30,6)))
                   #mut = c(6:8))

# for square matrix filled with 0.1s
design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0 # set diagonal to 0s

design

# run the feature selection tuning
tune.cll = tune.block.splsda(X = X, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 2, nrepeat = 1,
                              dist = "centroids.dist")

##### Results 


list.keepX = tune.cll$choice.keepX # set the optimal values of features to retain
list.keepX



####

final.diablo.model = block.splsda(X = X, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)

