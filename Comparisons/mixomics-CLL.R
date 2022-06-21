#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("mixOmics")

outdir2<-'Comparisons/plots/mixomics/'
#browseVignettes("mixOmics")

library(mixOmics)


data(liver.toxicity)
CLL_data

X <- liver.toxicity$gene


MyResult.pca <- pca(X)     # 1 Run the method
plotIndiv(MyResult.pca)    # 2 Plot the samples

plotVar(MyResult.pca)      # 3 Plot the variables

plotIndiv(MyResult.pca, group = liver.toxicity$treatment$Dose.Group, 
          legend = TRUE)



##### N integraiton 

library(mixOmics)
data(breast.TCGA)
# extract training data and name each data frame
X <- list(mRNA = breast.TCGA$data.train$mrna, 
          miRNA = breast.TCGA$data.train$mirna, 
          protein = breast.TCGA$data.train$protein)
Y <- breast.TCGA$data.train$subtype
summary(Y)

X <- list(mRNA = t(CLL_data$mRNA), 
            drug = t(CLL_data$Drugs))
Y <- CLL_metadata$IGHV

# FILTER MISSING VALUES 

rownames(CLL_data$mRNA)
tokeep<-which(!is.na(CLL_metadata$IGHV))
X <- list(mRNA = t(CLL_data$mRNA[,tokeep]), 
          drug = t(CLL_data$Drugs[,tokeep]), 
          meth=t(CLL_data$Methylation[,tokeep]),
          mut=t(CLL_data$Mutations[,tokeep])
          )
          
Y<-as.factor(CLL_metadata$IGHV[tokeep])


#list.keepX <- list(mRNA = c(16, 17), miRNA = c(18,5), protein = c(5, 5))
list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), drug=c(15,15))

list.keepX <- list(mRNA = c(16, 17), drug = c(18,5), meth=c(16,16), mut=c(15,15))

total_comps=15
MyResult.diablo <- block.splsda(X, Y, keepX=list.keepX,ncomp = total_comps)
plotIndiv(MyResult.diablo)



plotVar(MyResult.diablo, var.names = c(TRUE, TRUE, TRUE),
        legend=TRUE, pch=c(16,16,16))

selectVar(MyResult.diablo, block = 'mRNA', comp = 1)$mRNA$name



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
p<-circosPlot(MyResult.diablo, cutoff=0.7)
save(file= 'circos.png' )


##### cimDiablo
# Represent the multi-omics molecular signature expression for each sample 



#####

# minimal example with margins improved:
# cimDiablo(MyResult.diablo, margin=c(8,20))

# extended example:
cimDiablo(MyResult.diablo, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8,20), legend.position = "right")



###### Plot loadings: visualizes loading weights of each selected variables on each component
# And each dataset 
#plotLoadings(MyResult.diablo, contrib = "max")
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

cutoff=0.81
network(MyResult.diablo, blocks = c(1,3,4),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        lwd.edge = 2,
        cutoff = cutoff, 
        cex.node.name	=0.8,
        
        save = 'jpeg', name.save = paste0(outdir2, 'DIABLOnetwork_', cutoff))

