
### Load the rnas and mirnas

# 1. 
#BiocManager::install('OmnipathR')
#BiocManager::install('OmnipathR')

library('OmnipathR')

library(OmnipathR)
library(tidyr)
library(dnet)

#library(gprofiler2)
# interactions - proteins

mofa_cluster_id<-1
# TODO: fix to choose the visit ro plor.... 
# and the backbone visit 
VISIT_COMP<-'V08'
as.numeric(cell_corr)
cell_corr=TRUE
outdir
load_de_by_visit_and_cluster<-function(VISIT_COMP , mofa_cluster_id, cell_corr_state = FALSE, prefix='rnas_' ){
     #'
     #' @param VISIT_COMP
     #' @param mofa_cluster_id   
    #' @param cell_corr   

    rnas_visit<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c', as.numeric(cell_corr_state), '/' , VISIT_COMP ,  '/', prefix,  'de_cluster_', mofa_cluster_id, '.csv'))

    if (prefix == 'mirnas_'){
      rnas_visit$GENE_SYMBOL<-rnas_visit$X 
    }

    return(rnas_visit )

}

rnas_V08<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id,   cell_corr_state = cell_corr )
rnas_sig_V08<-rnas_V08%>% dplyr::filter(mofa_sign == 'Significant')

rnas_V08_list<-lapply(c(1,2,3), load_de_by_visit_and_cluster, VISIT_COMP='V08', cell_corr=cell_corr, prefix= 'rnas_') # holds all clusters 
rnas_sig_V08_list<-lapply(rnas_V08_list, function(df){ return(df%>% dplyr::filter(mofa_sign == 'Significant'))})


rnas_V06_list<-lapply(c(1,2,3), load_de_by_visit_and_cluster, VISIT_COMP='V06', cell_corr=cell_corr, prefix= 'rnas_') # holds all clusters 
rnas_sig_V06_list<-lapply(rnas_V06_list, function(df){ return(df%>% dplyr::filter(mofa_sign == 'Significant'))})

# TODO:  union of all visits 




## load mirs
mirnas_V08_list<-lapply(c(1,2,3), load_de_by_visit_and_cluster, VISIT_COMP=VISIT_COMP, cell_corr=cell_corr, prefix= 'mirnas_')
mirnas_sig_V08_list<-lapply(mirnas_V08_list, function(df){  return(df%>% dplyr::filter(mofa_sign == 'Significant'))})


mirnas_V06_list<-lapply(c(1,2,3), load_de_by_visit_and_cluster, VISIT_COMP='V06', cell_corr=cell_corr, prefix= 'mirnas_')
mirnas_sig_V06_list<-lapply(mirnas_V06_list, function(df){  return(df%>% dplyr::filter(mofa_sign == 'Significant'))})

mirnas_sig_V06_list





#### CHOOSE VISIT and cluster HERE
mirnas_visit<-mirnas_V08_list[[mofa_cluster_id]] # chosen visit to analyse
mirnas_sig_visit<-mirnas_sig_V08_list[[mofa_cluster_id]] # chosen visit to analyse

### define the rnas fcs that will be used
 
rnas_visit <- rnas_V08_list[[mofa_cluster_id]]
rnas_sig_visit<-rnas_sig_V08_list[[mofa_cluster_id]]


##### filter mofa top

rnas_top<-rnas_visit %>%arrange(padj) %>%  dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame
mirnas_top<-mirnas_visit %>%arrange(padj) %>% dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame


sel_factor=9; top_fr=0.05
top_genes_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'RNA',  factors=sel_factor, top_fr = top_fr))
top_genes_factor8<-get_symbols_vector(top_genes_factor8)
top_fr_mirs=top_fr
top_mirnas_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'miRNA',  factors=sel_factor, top_fr = top_fr_mirs))


## 1. Colour if in TOP 1% of factor
## 2. Colour if in specific pathways 
## 
dim(rnas_top); dim(mirnas_top)
top_g<-10000
top_m<-200

top_rnas<-rnas_top[1:top_g,]
top_mirnas<-mirnas_top[1:top_m,]
top_mirnas$GENE_SYMBOL<-top_mirnas$X






source(paste0(script_dir, 'ppmi/network_utils.R'))
############# MY ADDITION 
# 1. DE MIRS 
# 2. DE GENES 
mirnas_sig_visit

mirnas_sig_factor<-mirnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
rnas_sig_factor<-rnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)


## compare clust 1,2

mirnas_sig_factor_V08_list<-lapply(mirnas_sig_V08_list, function(df){
  return( df %>%  dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8))
})

rnas_sig_factor_V08_list<-lapply(rnas_sig_V08_list, function(df){
  return( df %>%  dplyr::filter(GENE_SYMBOL %in% top_genes_factor8))
})

## input: 
## de_rnas: union of all
## de_mirnas: union of all clusters 
rnas_sig_factor_V08_list
rnas_sig_factor_all_clusts<-Reduce(function(x, y) merge(x, y, all=TRUE), rnas_sig_factor_V08_list)
mirnas_sig_factor_all_clusts<-Reduce(function(x, y) merge(x, y, all=TRUE), mirnas_sig_factor_V08_list)

mirnas_sig_factor_all_clusts$GENE_SYMBOL
g_extended_both_clusters<-create_regulatory_net_backbone(rnas_sig_factor_all_clusts, mirnas_sig_factor_all_clusts)

g_extended_both_clusters
#g_extended<-create_regulatory_net_backbone(rnas_sig_factor, mirnas_sig_factor)

### 
g_extended_both_clusters<-create_regulatory_net_backbone_mirtar(rnas_sig_factor_all_clusts, mirnas_sig_factor_all_clusts)


### Now plot! 
g_extended_both_clusters
OPI_g_union<-g_extended_both_clusters
# TODO: give the layout 

set.seed(123) 
# color by logFC 

# backbone should be the union 
#g_fc_V08<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08,  de_mirnas=mirnas_V08_list[[3]])




# TODO: apply take mis and genes 

mofa_cluster_id

rnas_V08_list[[2]]

mirnas_V08_list[[1]]
# Unified plot: 
g_fc_V08_list<-sapply(c(1,2,3), function(mofa_cluster_id){
  g_fc_V08_cl<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08_list[[mofa_cluster_id]],  de_mirnas=mirnas_V08_list[[mofa_cluster_id]])
  return(g_fc_V08_cl)

})

g_fc_V06_list<-sapply(c(1,2,3), function(mofa_cluster_id){
  g_fc_V06_cl<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V06_list[[mofa_cluster_id]],  de_mirnas=mirnas_V06_list[[mofa_cluster_id]])
  return(g_fc_V06_cl)

})
g_fc_V06_list

OPI_g_union




## igraph plotting 

### Plot with specified coords 
# Get the coordinates of the Nodes of the backbone that includes all clusters 

Coords <- layout_with_fr(OPI_g_union)# %>% 
  #as_tibble %>%
  #  bind_cols(data_frame(names = names(V(g_fc))))



## Plot cluster 3 as specified by backbone coords 

dir.create(paste0(outdir, '/networks/'))
mofa_cluster_id=1
set.seed(123) 

g_fc_plot<-g_fc_V08_list[[mofa_cluster_id]];VISIT_PLOT = 'V08'

g_fc_plot<-g_fc_V06_list[[mofa_cluster_id]];  VISIT_PLOT = 'V06'

## ADD border if significant 
# create a vertor of border colours conditional on node type

V(g_fc_V06_list[[mofa_cluster_id]])$FC


V(g_fc_V08_list[[mofa_cluster_id]])$FC

round(cbind(V(g_fc_V08_list[[mofa_cluster_id]])$FC,V(g_fc_V06_list[[mofa_cluster_id]])$FC), digits=2)

sig_names<-V(g_fc_plot)$name[V(g_fc_plot)$significant]
bd <- ifelse(V(g_fc_plot)$significant, "#2fc729", NA) 
V(g_fc_plot)$size<-ifelse(!is.na(V(g_fc_plot)$FC), log2(abs(V(g_fc_plot)$FC)+1)*17, 0.1*20)
g_fc_plot$layout<-Coords

# IF Using omnipath

#V(g_fc_plot)$size<-V(g_fc_plot)$size*
plot_settings_omnipath<-list(
  label.cex=2,
  vertex.size.factor=1.5 
)
plot_settings_mirtar= list(
  label.cex=1,
  vertex.size.factor=1 
)


plot_settings=plot_settings_mirtar

V(g_fc_plot)$size<-V(g_fc_plot)$size*plot_settings$vertex.size.factor

 ### IGRAPH TO SAVE 

net_name=paste0('mirs_genes_',VISIT_PLOT,'_', mofa_cluster_id, '_f',sel_factor,top_fr )
V(g_fc_plot)$frame.color<-bd
V(g_fc_plot)$label.cex=2 # increase label size? 
V(g_fc_plot)$label.cex=plot_settings$label.cex # increase label size? 


V(g_fc_plot)
plot(g_fc_plot, 
vertex.color=(adjustcolor(V(g_fc_plot)$color, alpha.f = 0.8 )), 
vertex.size=V(g_fc_plot)$size,
vertex.frame.color=bd,
vertex.frame.width=3,
 vertex.label.dist = 2, 
  layout = Coords)
title(paste('Visit: ', VISIT_COMP, ', Cluster: ', mofa_cluster_id, ', Factor: ',
 sel_factor, top_fr  ),cex.main=3,col.main="green")





#### VISNWT TO SAVE


#V(g_fc)$shape<- ifelse(grepl("miR|hsa-let",igraph::V(g)$name), "vrectangle", "circle")
visnet <- toVisNetworkData(g_fc_plot)
visnet$nodes$abs_FC<- abs(visnet$nodes$FC)
 
visSave(visIgraph(g_fc_plot), file = paste0(outdir, '/networks/',  net_name, '.html'))

visIgraph(g_fc_plot, layout=Coords) 








write.csv(V(g_fc_plot)$name,paste0(outdir, '/genes_network.csv'))


 visualize_net(visnet)
visnet<-toVisNetworkData(g_fc_plot, idToLabel = TRUE)
bd <- ifelse(V(g_fc_plot)$significant, "yellow", NA) 

visNetwork(visnet$nodes, visnet$edges)%>% 
  visNodes(color = list(background = "blue", 
  border = "black"), borderWidth = 3)
























































