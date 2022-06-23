

## Choose here the model to run, supervised or not 

MyResult.diablo<-final.diablo.model
total_comps=ncomp
dims_x<-length(names(X))
params_str<-paste0('sup_',paste(unlist(names(X)), collapse='_'), '_', length(names(X)), '_', total_comps)

png(paste0 (outdir2, 'feature_plot', params_str,'.png'))
plotVar(MyResult.diablo, var.names = c(TRUE, TRUE, TRUE),
        legend=TRUE, pch=c(16,16,16))

dev.off()
# print the first component 
selectVar(MyResult.diablo, block = 'mRNA', comp = 1)$mRNA$name
selectVar(MyResult.diablo, block = 'meth', comp = 1)$meth$name


##### Customize sample plots 
png(paste0(outdir2,'factor_space_', params_str,'.png'))
plotIndiv(MyResult.diablo, 
          ind.names = FALSE, 
          legend=TRUE, cex=c(1,2),
          title = 'CLL with DIABLO')

dev.off()

# Only show the names of the proteins 
plotVar(MyResult.diablo, var.names = rep(TRUE, dims_x),
        legend=TRUE, pch=rep(16,dims_x ))

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

##### GLoval OVERVIEW of the correlation structure at the componet level 
png(paste0(outdir2, 'correlation_structure', length(names(X)), '.png' ))
plotDiablo(MyResult.diablo, ncomp = 1)
dev.off()

png(paste0(outdir2, 'arrow', params_str, '.png' ))
plotArrow(MyResult.diablo, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
dev.off()

### Visualize correlations between variables 
p<-circosPlot(MyResult.diablo, cutoff=0.9)
save(file= 'circos.png' )


##### cimDiablo
# Represent the multi-omics molecular signature expression for each sample 



#####

# minimal example with margins improved:
cimDiablo(MyResult.diablo, margin=c(8,20))

# extended example:
color.blocks=c('darkorchid', 'brown1', 'lightgreen','yellow' )[1:dims_x]
png(paste0(outdir2,'cim_', params_str, '.png'))
cimDiablo(MyResult.diablo, color.blocks = color.blocks,
           margin=c(8,20), legend.position = "right" 
)
dev.off()




###### Plot loadings: visualizes loading weights of each selected variables on each component
# And each dataset \
png(paste0(outdir2,'loadings', params_str, '_', comp, '.png'))
plotLoadings(MyResult.diablo, contrib = "max")
dev.off()
comp=1
p<-plotLoadings(MyResult.diablo, comp = comp, contrib = "max")
show(p)
ggsave(paste0(out_dir,'top_weights', comp, '.png'), width = 4, height = 4, dpi=500 )
dev.off()




P<-plotLoadings(MyResult.diablo, comp = 1, contrib = "max")
ggsave('top_weights.png', width = 4, height = 4, dpi=500 )

comps=c(1:5)





v_set=c()

for (comp in 1:total_comps){ 
  top_genes_mixomics<-selectVar(MyResult.diablo, comp = comp,list.keepX=10)$mRNA$name[1:list.keepX$mRNA[1]]
  write.csv(top_genes_mixomics, paste0(outdir2, 'top_genes_', total_comps, '_', comp, '.csv'))
  
  top_meth_mixomics<-selectVar(MyResult.diablo, comp = comp,list.keepX=10)$meth$name[1:list.keepX$meth[1]]
  write.csv(top_meth_mixomics, paste0(outdir2, 'top_meth_', total_comps, '_', comp, '.csv'))
  
  
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
blocks=c(1,2,3,4)[1:dims_x]

#png(paste0(outdir2, 'diablo_plot',param_str, '.png'))
network(MyResult.diablo, blocks = blocks,
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        cutoff = 0.85, save = 'jpeg', name.save = 'DIABLOnetwork')

#dev.off()

cutoff=0.85
network(MyResult.diablo, blocks = blocks,
        color.node = c('darkorchid', 'brown1', 'lightgreen'), 
        lwd.edge = 2,
        cutoff = cutoff, 
        cex.node.name	=0.8,
        save = 'jpeg', name.save = paste0(outdir2, 'DIABLOnetwork_', as.character(params_str),'_', paste0(blocks, collapse='_'), cutoff))

