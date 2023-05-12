
#BiocManager::install(c("edgeR"))
#BiocManager::install(c("Glimma"))

library(mixOmics)
#### 
# Preprocessing 


#first remember the names
library(dplyr )
library(resample)
library(mixOmics)

## NEEDED FOR CPM function in gene preprocessing" 
#install.packages('gplot')
library(edgeR)
library(limma)
library(Glimma)


ncomp_g=5
ncomp_p=5
ncomp=5


X1_t<-t(assays(mofa_multi_rna_mir_complete)[['RNA']])
X2_t<-t(assays(mofa_multi_rna_mir_complete)[['miRNA']])
mofa_multi_rna_mir_complete

dim(X1_t)
dim(X2_t)
rownames(X1_t)==rownames(X2_t)

result.pls <- pls(X1_t, X2_t, ncomp = ncomp, mode='canonical')  # where ncomp is the number of dimensions/components to choose
perf.spls <- perf(result.pls, validation = "loo",
                 progressBar = FALSE, nrepeat = 10)








