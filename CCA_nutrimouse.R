
## install BiocManager if not installed if (!requireNamespace("BiocManager", quietly = TRUE))     install.packages("BiocManager") 
## install mixOmics 
BiocManager::install('mixOmics', force=TRUE)

library(mixOmics) # import the mixOmics library

set.seed(5249) # for reproducibility, remove for normal use

# Case Study of rCCA with Nutrimouse dataset
data(nutrimouse)
X <- nutrimouse$lipid # extract all lipid concentration variables
Y <- nutrimouse$gene # extract all gene expression variables

dim(X) # check the dimensions of the X dataf



## exploration of feature selection


imgCor(X, Y, sideColors = c("purple", "green")) 



# set grid search values for each regularisation parameter
grid1 <- seq(0.001, 0.2, length = 10) 
grid2 <- seq(0.001, 0.2, length = 10)

# optimise the regularisation parameter values
cv.tune.rcc.nutrimouse <- tune.rcc(X, Y, grid1 = grid1, grid2 = grid2, 
                                   validation = "loo") 



cv.tune.rcc.nutrimouse # examine the results of CV tuning


opt.l1 <- cv.tune.rcc.nutrimouse$opt.lambda1 # extract the optimal lambda values
opt.l2 <- cv.tune.rcc.nutrimouse$opt.lambda2

# formed optimised CV rCCA
CV.rcc.nutrimouse <- rcc(X, Y, method = "ridge", 
                         lambda1 = opt.l1, lambda2 = opt.l2) 



# run the rCCA method using shrinkage
shrink.rcc.nutrimouse <- rcc(X,Y, method = 'shrinkage') 
# examine the optimal lambda values after shrinkage 
shrink.rcc.nutrimouse$lambda 



# barplot of cross validation method rCCA canonical correlations
plot(CV.rcc.nutrimouse, type = "barplot", main = "Cross Validation") 

# barplot of shrinkage method rCCA canonical correlations
plot(shrink.rcc.nutrimouse, type = "barplot", main = "Shrinkage") 



# plot the projection of samples for CV rCCA data
plotIndiv(CV.rcc.nutrimouse, comp = 1:2, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) Nutrimouse, rCCA CV XY-space')

# plot the projection of samples for shrinkage rCCA data
plotIndiv(shrink.rcc.nutrimouse, comp = 1:2, 
          ind.names = nutrimouse$genotype,
          group = nutrimouse$diet, rep.space = "XY-variate", 
          legend = TRUE, title = '(b) Nutrimouse, rCCA shrinkage XY-space')




# plot the arrow plot of samples for CV rCCA data
plotArrow(CV.rcc.nutrimouse, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method')

# plot the arrow plot of samples for shrinkage rCCA data
plotArrow(shrink.rcc.nutrimouse, group = nutrimouse$diet, 
          col.per.group = color.mixo(1:5),
          title = '(b) Nutrimouse, shrinkage method')


