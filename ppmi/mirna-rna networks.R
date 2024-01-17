
library(igraph)
merged_targets




### TODO: RUN AFTER MIRNA SEQ ENRICHMENT WITH CORRELATION 

### Filter also by Differentially expressed genes only??? 
deseq2ResDF = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)
process_mirnas=TRUE; source(paste0(script_dir, '/config.R'))
mirs


gene_list


outdir_s_mirnas<-outdir_s

process_mirnas=FALSE; source(paste0(script_dir, '/config.R'))

de_rnas = read.csv(paste0(outdir_s, '/significant0.05_0.1.csv'), row.names = 1)

de_rnas_all = read.csv(paste0(outdir_s, '/results_df.csv'), row.names = 1)
de_rnas_sig<-de_rnas_all[de_rnas_all$padj<0.05 & de_rnas_all$log2FoldChange>0.05 ,]
rownames(de_rnas_sig)
de_ids<-gsub('\\..*', '',rownames(de_rnas_sig))


merged_targets_onlyde<-merged_targets[merged_targets$target_ensembl %in% rownames(de_rnas),]
merged_targets_onlyde<-merged_targets[merged_targets$target_ensembl %in% de_ids,]




mieaa_targets_de_rnas$Subcategory
dcast(mieaa_targets_de_rnas$miRNAs.precursors,  )


merged_targets_el<-merged_targets_onlyde[, c('mature_mirna_id', 'symbol', 'cor')]

# definitions 
# 1. de mirs
# 2. de genes 
de_mirs<-mirnas_V08_cl1
de_mirs<-mirnas_V08_cl1
mirna_targets_el<-merged_targets_el
de_rnas_all<-rnas_V08_cl1$GENE_SYMBOL


############

g<-graph_from_edgelist(as.matrix(mirna_targets_el[,c(1,2)]), directed = TRUE)
g<-graph_from_edgelist(as.matrix(mirna_targets), directed = TRUE)

#E(g)$weight <- merged_targets_onlyde[,'cor']
g$layout <- layout_with_kk

p<-plot.igraph(g )






ggsave(paste0(outdir_s_mirnas, '/mirna_rna_net.svg'), plot=last_plot())




##### 
mst<-g

mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1
V(mst) %in% rownames(de_mirs)
mir_ids<-!is.na(match(names(V(mst)), rownames(de_mirs)))
new=c()
de_mirs[ match(names(V(mst)), rownames(de_mirs) ),]$log2FoldChange
new=de_mirs[ match(names(V(mst)), rownames(de_mirs) ),]$log2FoldChange
de_mirs[rownames(de_mirs)%in%c('hsa-miR-26b-5p'),]
id_rnas<-match(names(V(mst)), de_rnas_all$SYMBOL)
new

rnas_lfc<-de_rnas_all[ id_rnas ,]$log2FoldChange
rnas_lfc<-rnas_lfc[!is.na(rnas_lfc)]
rnas_lfc
new[!is.na(id_rnas)]<-rnas_lfc
V(g)$name
vSizes <-new
V(mst)$node
deseq2ResDF$log2FoldChange



par(mfrow=c(1,2))
plot(
  mst.clustering, mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  #vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.7,
 # edge.width=edgeweights,
  edge.arrow.mode=0,
  main="RNA-MIRNA anticorelated targets")

plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.7,
  #edge.width=edgeweights,
  edge.arrow.mode=0,
 main="RNA-MIRNA anticorelated targets"
)

dev.off()
plot(
 # mst.clustering,
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  #vertex.size=abs(vSizes)*25,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  
  asp=FALSE,
  vertex.label.cex=1,
  # edge.width=edgeweights,
  edge.arrow.mode=0,
  main="RNA-MIRNA anticorelated targets")


ggsave(paste0(outdir_s_mirnas, '/mirna_rna_net.svg'), plot=last_plot())



library(visNetwork)
library(htmlwidgets)
#install.packages('visNetwork')
#install.packages('htmlwidgets')


saveWidget(visIgraph(g), file = paste0(outdir_s_mirnas, '/mirna_rna_net.html'))











