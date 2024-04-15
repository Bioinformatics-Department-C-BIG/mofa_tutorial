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


list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), mut=c(15,15))

# without mutations
list.keepX <- list(mRNA = c(15, 17), drug = c(18,5), meth=c(15,5))

total_comps=5

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





plotVar(MyResult.diablo, var.names = c(TRUE, TRUE, TRUE),
        legend=TRUE, pch=c(16,16,16))


# print the first component 
selectVar(MyResult.diablo, block = 'mRNA', comp = 1)$mRNA$name
selectVar(MyResult.diablo, block = 'meth', comp = 1)$meth$name


##### Customize sample plots 

plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2),
          title = 'CLL with DIABLO')




dev.off()
# Only show the names of the proteins 
plotVar(MyResult.diablo, var.names = c(TRUE, TRUE),
        legend=TRUE, pch=c(16,16))

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
cimDiablo(MyResult.diablo, margin=c(8,20))

# extended example:
png(paste0(outdir2,'cim.png'))
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen','yellow' ),
          comp = 1, margin=c(8,20), legend.position = "right", 
          )
dev.off()


###### Plot loadings: visualizes loading weights of each selected variables on each component
# And each dataset 
#plotLoadings(MyResult.diablo, contrib = "max")
comp=1
p<-plotLoadings(MyResult.diablo, comp = comp, contrib = "max")
ggsave(paste0(out_dir,'top_weights', comp, '.png'), width = 4, height = 4, dpi=500 )
dev.off()




P<-plotLoadings(MyResult.diablo, comp = 1, contrib = "max")
ggsave('top_weights.png', width = 4, height = 4, dpi=500 )

comps=c(1:5)





v_set=c()

for (comp in 1:15){ 
  top_genes_mixomics<-map_to_gene(selectVar(MyResult.diablo, comp = comp,list.keepX=20)$Mutations$name)$SYMBOL[1:15]
  write.csv(top_genes_mixomics, paste0(outdir2, 'top_genes_', total_comps, '_', comp, '.csv'))
  
  # concat all 
  v_set <- c(v_set, top_genes_mixomics)
  print(length(v_set))
}
write.csv(v_set, paste0(outdir2, 'all_mutations_top_', total_comps, '.csv'))

#write.csv(v_set, paste0(outdir2, 'all_genes_top_', total_comps, '.csv'))

v_set
outdir2

ens_ids2
?network

network(MyResult.diablo, blocks = c(1,2,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.85, save = 'jpeg', name.save = 'DIABLOnetwork')

dev.off()

cutoff=0.80
network(MyResult.diablo, blocks = c(1,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        lwd.edge = 2,
        cutoff = cutoff, 
        cex.node.name	=0.8,
        save = 'jpeg', name.save = paste0(outdir2, 'DIABLOnetwork_', cutoff))

