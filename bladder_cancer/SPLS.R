
#### 
# Preprocessing 


#first remember the names
library(dplyr )
install.packages('resample')
library(resample)
library(mixOmics)

transpose_matrix<- function(df.aree){
    n <- df.aree$Symbol
    # transpose all but the first column (name)
    df.aree <- as.data.frame(t(df.aree[,-1]))
    colnames(df.aree) <- n
    #df.aree$myfactor <- factor(row.names(df.aree))
  return(df.aree)
    
  }


preprocess_raw_data<-function(df){
      df<- as.data.frame(apply(df, 2, function(x) as.numeric(x)) )
      ind <- apply(df, 2, function(x) sum(x, na.rm = TRUE)==0) 
      df<-df[,!ind]
      return(df)

}

####################################
#### Preliminary analysis with PCA


data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
head(cbind(rownames(X), rownames(Y)))

pca.gene <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)



plot(pca.gene)

dir='H:/My Drive/PHD 2020/Projects/Bladder cancer/'
X1_raw<-read.csv(file = paste0(dir,'RNAseq_BladderCancer.csv' ))
X2_raw<-read.csv(file = paste0(dir,'Proteomics_BladderCancer.csv' ))
Y_raw<-read.csv(file = paste0(dir,'pheno_BladderCancer.csv' ), nrows = 16)

X1_t<-transpose_matrix(X1_raw)
X1_t<-X1_t[,1:27000]
df<-X1_t
df<-as.data.frame(apply(df, 2, function(x) as.numeric(x)))
ind <- apply(df, 2, function(x) sum(x, na.rm = TRUE)==0) 
#remove genes with zero variance
df<-df[,!ind]
X1_t<-df
#X1_t<-preprocess_raw_data(X1_t)


Y_raw$Subtype<-as.factor(Y_raw$Subtype)
################
#### Preprocessing -transpose


pca.gene <- pca(X1_t, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)




plotIndiv(pca.gene, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')


spca.result <- spca(X1_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 5, 15))

selectVar(spca.result, comp = 1)$value

####################################################################
##### SPLS DA

library(mixOmics)
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic


pca.clinical <- pca(Y, ncomp = 10, center = TRUE, scale = TRUE)
head(selectVar(pca.gene, comp=1)$value)

plotIndiv(pca.gene, comp = c(1, 2), group = liver.toxicity$treatment[, 4],
          ind.names = liver.toxicity$treatment[, 3],
          legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')


#PLS
result <- pls(X, Y, ncomp = 3)  # where ncomp is the number of dimensions/components to choose
tune.pls <- perf(result, validation = 'loo', criterion = 'all', progressBar = FALSE)


# SPLS
ncomp = 10
result.spls <- spls(X, Y, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
                  criterion = 'all', progressBar = FALSE)






library(mixOmics)
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene


nutrimouse.shrink <- rcc(X, Y, ncomp = 3, method = 'shrinkage')


plot(nutrimouse.shrink, type = "barplot")


## Estimation of penalisation parameters (CV method)
grid1 <- seq(0, 0.2, length = 5) 
grid2 <- seq(0.0001, 0.2, length = 5)

cv <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, validation = "loo")


result <- rcc(X, Y, ncomp = 3, lambda1 = cv$opt.lambda1, lambda2 = cv$opt.lambda2)
result$cor


