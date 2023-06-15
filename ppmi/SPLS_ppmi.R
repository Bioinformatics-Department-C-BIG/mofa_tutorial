
#BiocManager::install(c("edgeR"))
#BiocManager::install(c("mixOmics"))

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




########### DIABLO



data = list(RNA = X1_t, 
            miRNA = X2_t )

Y<- mofa_multi_rna_mir_complete$COHORT_DEFINITION



design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0

design 

# tune components
sgccda.res = block.splsda(X = data, Y = Y, ncomp = 7, 
                          design = design)

set.seed(123) # for reproducibility, only when the `cpus' argument is not used
# this code takes a couple of min to run

## This chooses the components 
perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = 4, 
                   nrepeat = 4, progressBar = TRUE)

#perf.diablo  # lists the different outputs
plot(perf.diablo) 

perf.diablo$choice.ncomp$WeightedVote
ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]

ncomp

## This chooses the number of features 
test.keepX <- list(RNA = c(6:9, seq(10, 25, 10)),
                   miRNA = c(6:9, seq(10, 20, 5)))

tune.diablo <- tune.block.splsda(data, Y, ncomp = ncomp, 
                                      test.keepX = test.keepX, design = design,
                                      validation = 'Mfold', folds = 3, nrepeat = 1, 
                                      BPPARAM = BiocParallel::SnowParam(workers = 3),
                                      dist = "centroids.dist")

list.keepX <- tune.diablo.tcga$choice.keepX



### SHOW most important variables 
sgccda.res$design
selectVar(sgccda.res, block = 'RNA', comp = 1)$RNA$name 




plotDiablo(perf.diablo, ncomp = 1)
plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')




