
#### 
# Preprocessing 


#first remember the names
library(dplyr )
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

dir='/Users/efiathieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'
dir='E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Projects/Bladder cancer/'

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

X2_t<-transpose_matrix(X2_raw)
df<-X2_t
df<-as.data.frame(apply(df, 2, function(x) as.numeric(x)))
ind <- apply(df, 2, function(x) sum(x, na.rm = TRUE)==0) 
#remove genes with zero variance
df<-df[,!ind]
X2_t<-df



Y_raw$Subtype<-as.factor(Y_raw$Subtype)
################
#### Preprocessing -transpose


pca.gene <- pca(X1_t, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.gene)
plotIndiv(pca.gene, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')


spca.result <- spca(X1_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 10,5))

selectVar(spca.result, comp = 1)$value




pca.proteomics <- pca(X2_t, ncomp = 10, center = TRUE, scale = TRUE)
plot(pca.proteomics)
plotIndiv(pca.proteomics, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')


spca.result <- spca(X2_t, ncomp = 3, center = TRUE, scale = TRUE, 
                    keepX = c(10, 10, 5))

selectVar(spca.result, comp = 2)$value





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
result.spls <- spls(X1_t, X2_t, ncomp = ncomp, keepX = c(rep(10, ncomp)), mode = 'regression')
tune.spls <- perf(result.spls, validation = 'Mfold', folds = 10,
                  criterion = 'all', progressBar = FALSE)


plotIndiv(pca.proteomics, comp = c(1, 2), group = Y_raw$Subtype,
          legend = TRUE, title = 'Liver gene, PCA comp 1 - 2')


plotIndiv(result.spls, ind.names = FALSE,
          rep.space = "XY-variate", # plot in Y-variate subspace
          group =  Y_raw$Subtype, # colour by time group
          legend = TRUE)


# Enrichment analysis

reticulate::use_python(python = "C:/Users/athienitie/Anaconda3/python.exe")


library(mixOmics)
data(breast.TCGA)
#BiocManager::install("biomaRt")

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)

#### Enrichment analysis 
utils::data()



# GSEA on positive weights, with default options
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "positive"
)



#DIABLO
data = list(mRNA = X1_t, 
            proteomics = X2_t )

Y<-Y_raw$Subtype


design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

design 

# tune components
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 5, 
                         design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, nrepeat = 10)

#perf.diablo  # lists the different outputs
plot(perf.diablo) 


sgccda.res$design
selectVar(sgccda.res, block = 'mRNA', comp = 1)$mRNA$name 


plotDiablo(sgccda.res, ncomp = 1)
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')



plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

dev.off()  
plotVar(sgccda.res, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'brown1' ))


circosPlot(sgccda.res, cutoff = 0.7, line = TRUE, 
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3"), size.labels = 1.5)
