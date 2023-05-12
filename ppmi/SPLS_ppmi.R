
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
ncomp_p=4
ncomp=2


X1_t<-t(assays(mofa_multi_rna_mir_complete)[['RNA']])
X2_t<-t(assays(mofa_multi_rna_mir_complete)[['miRNA']])
mofa_multi_rna_mir_complete

dim(X1_t)
dim(X2_t)
rownames(X1_t)==rownames(X2_t)

result.pls <- pls(X1_t, X2_t, ncomp = ncomp, mode='regression')  # where ncomp is the number of dimensions/components to choose

perf.spls <- perf(result.pls, validation = "Mfold",
                 progressBar = TRUE,
                 folds=2, nrepeat = 10)


#png(paste0(output,'/PLS_Q2',param_str,'.png'))
plot(perf.spls, criterion = 'Q2.total')
#abline(h = 0.0975)
#dev.off()



#### TUNING
# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))
# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10) 


tune.spls.ppmi <- tune.spls(X1_t, X2_t, ncomp = 2,
                             test.keepX = list.keepX,
                             test.keepY = list.keepY,
                             nrepeat = 3, folds = 2, # use 10 folds
                             mode = 'regression', measure = 'cor',
                            progressBar = TRUE) 
plot(tune.spls.ppmi)    


optimal.keepX <- tune.spls.ppmi$choice.keepX 
optimal.keepY <- tune.spls.ppmi$choice.keepY
optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components


# use all tuned values from above
final.spls.liver <- spls(X1_t,X2_t , ncomp = optimal.ncomp, 
                         keepX = optimal.keepX,
                         keepY = optimal.keepY,
                         mode = "regression") # explanitory approach being used, 
# hence use regression mode

plotIndiv(final.spls.liver, ind.names = FALSE, 
          rep.space = "X-variate", # plot in X-variate subspace
          group = mofa_multi_rna_mir_complete$COHORT_DEFINITION, # colour by time group
          #pch = as.factor(liver.toxicity$treatment$Dose.Group), 
          col.per.group = color.mixo(1:2), 
          legend = TRUE, legend.title = 'Time', legend.title.pch = 'Dose')






