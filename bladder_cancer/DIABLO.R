library('mixOmics')

data(breast.TCGA) # load in the data


# set a list of all the X dataframes
data = list(mRNA = X1_t,
            proteomics = X2_t)

lapply(data, dim) # check their dimensions
Y<-Y_raw$EORTC.risk

summary(Y)

list.keepX = c(20, 20) # select arbitrary values of features to keep
list.keepY = c(25, 2)

# generate three pairwise PLS models
pls1 <- spls(data[["mRNA"]], data[["proteomics"]], 
             keepX = list.keepX, keepY = list.keepY) 
dev.off()


# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) miRNA vs mRNA", 
        legend = c("mRNA", "proteomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# correlation btw the first components
cor(pls1$variates$X, pls1$variates$Y) 


design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

design

# form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) 



### tune the number of components
# run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 4, nrepeat = 10) 

plot(perf.diablo) # plot output of tuning



# set grid of values for each component to test
test.keepX = list (mRNA = c(6:7, seq(20,30,5)), 
                   proteomics = c(6:7,  seq(20,30,5)))

#run the feature selection tuning
tune.TCGA = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                       test.keepX = test.keepX, design = design,
                       validation = 'Mfold', folds = 10, nrepeat = 1,
                       dist = "centroids.dist")



#####
#list.keepX = list(mRNA=10,proteomics=10)
# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)

# the features selected to form the first component
head(selectVar(final.diablo.model, block = 'mRNA', comp = 2)$mRNA$name )
head(selectVar(final.diablo.model, block = 'proteomics', comp = 1)$proteomics$name )



##### PLOTS

plotDiablo(final.diablo.model, ncomp = 1)

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots')

plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')


plotVar(final.diablo.model, cutoff=0.8, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'brown1'))


circosPlot(final.diablo.model, cutoff = 0.98, line = TRUE,
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

