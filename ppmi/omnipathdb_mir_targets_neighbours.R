
### Load the rnas and mirnas

# 1. 

#BiocManager::install('OmnipathR')
library('OmnipathR')

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)
# interactions - proteins

mofa_cluster_id<-2
VISIT_COMP<-'V08'
as.numeric(cell_corr)
cell_corr=TRUE
load_de_by_visit_and_cluster<-function(VISIT_COMP , mofa_cluster_id, cell_corr_state = FALSE, prefix='rnas_' ){
     #'
     #' @param VISIT_COMP
     #' @param mofa_cluster_id   
    #' @param cell_corr   

    rnas_visit<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c', as.numeric(cell_corr_state), '/' , VISIT_COMP ,  '/', prefix,  'de_cluster_', mofa_cluster_id, '.csv'))

    return(rnas_visit )

}

rnas_V08<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id,   cell_corr_state = cell_corr )
rnas_sig_V08<-rnas_V08%>% dplyr::filter(mofa_sign == 'Significant')

rnas_V08_cl1<-load_de_by_visit_and_cluster(VISIT_COMP, 1, cell_corr_state = cell_corr )
rnas_V08_cl2<-load_de_by_visit_and_cluster(VISIT_COMP, 2 , cell_corr_state = cell_corr)
rnas_V08_cl3<-load_de_by_visit_and_cluster(VISIT_COMP, 3 , cell_corr_state = cell_corr )




rnas_sig_V08_cl2<-load_de_by_visit_and_cluster(VISIT_COMP, 2,  , cell_corr_state = cell_corr )%>% dplyr::filter(mofa_sign == 'Significant')
rnas_sig_V08_cl1<-load_de_by_visit_and_cluster(VISIT_COMP, 1 , cell_corr_state = cell_corr )%>% dplyr::filter(mofa_sign == 'Significant')
rnas_sig_V08_cl3<-load_de_by_visit_and_cluster(VISIT_COMP, 3 , cell_corr_state = cell_corr )%>% dplyr::filter(mofa_sign == 'Significant')


rnas_V06<-load_de_by_visit_and_cluster(VISIT_COMP='V06', mofa_cluster_id , cell_corr_state = cell_corr )
rnas_V04<-load_de_by_visit_and_cluster(VISIT_COMP='V04', mofa_cluster_id , cell_corr_state = cell_corr )
rnas_BL<-load_de_by_visit_and_cluster(VISIT_COMP='BL', mofa_cluster_id , cell_corr_state = cell_corr )

rnas_sig_V06<-rnas_V06%>% dplyr::filter(mofa_sign == 'Significant')



# TODO: union of all clusters? 
#

### define the rnas fcs that will be used 
rnas_visit <- rnas_V08
rnas_sig_visit<-rnas_sig_V08



## load mirs
mirnas_V08<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id = mofa_cluster_id,  cell_corr = cell_corr , prefix='mirnas_'); mirnas_V08$GENE_SYMBOL<-mirnas_V08$X
mirnas_visit<-mirnas_V08
mirnas_V08_cl1<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id = 1 , cell_corr = cell_corr, prefix='mirnas_'); mirnas_V08_cl1$GENE_SYMBOL<-mirnas_V08_cl1$X
mirnas_V08_cl2<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id = 2 , cell_corr = cell_corr, prefix='mirnas_'); mirnas_V08_cl2$GENE_SYMBOL<-mirnas_V08_cl2$X
mirnas_V08_cl3<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id = 3 , cell_corr = cell_corr, prefix='mirnas_'); mirnas_V08_cl3$GENE_SYMBOL<-mirnas_V08_cl3$X


mirnas_V06<-load_de_by_visit_and_cluster('V06', mofa_cluster_id = mofa_cluster_id , cell_corr = cell_corr, prefix='mirnas_'); mirnas_V06$GENE_SYMBOL<-mirnas_V06$X
mirnas_V04<-load_de_by_visit_and_cluster('V04', mofa_cluster_id = mofa_cluster_id , cell_corr = cell_corr, prefix='mirnas_'); mirnas_V04$GENE_SYMBOL<-mirnas_V04$X
mirnas_BL<-load_de_by_visit_and_cluster('BL', mofa_cluster_id = mofa_cluster_id , cell_corr = cell_corr, prefix='mirnas_'); mirnas_BL$GENE_SYMBOL<-mirnas_BL$X


mirnas_sig_V08_cl1<-mirnas_V08_cl1%>% dplyr::filter(mofa_sign == 'Significant')
mirnas_sig_visit<-mirnas_sig_V08_cl1
mirnas_sig_V08_cl2<-mirnas_V08_cl2%>% dplyr::filter(mofa_sign == 'Significant')
mirnas_sig_V08_cl3<-mirnas_V08_cl3%>% dplyr::filter(mofa_sign == 'Significant')


rnas_top<-rnas_visit %>%arrange(padj) %>%  dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame
mirnas_top<-mirnas_visit %>%arrange(padj) %>% dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame

sel_factor=8; top_fr=0.1
top_genes_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'RNA',  factors=sel_factor, top_fr = top_fr))
top_genes_factor8<-get_symbols_vector(top_genes_factor8)

top_fr_mirs=top_fr
top_mirnas_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'miRNA',  factors=sel_factor, top_fr = top_fr_mirs))


length(top_genes_factor8)


## 1. Colour if in TOP 1% of factor
## 2. Colour if in specific pathways 
## 
top_genes_factor8

dim(rnas_top); dim(mirnas_top)
top_g<-10000
top_m<-200

top_rnas<-rnas_top[1:top_g,]
top_mirnas<-mirnas_top[1:top_m,]

top_mirnas$GENE_SYMBOL<-top_mirnas$X
top_mirnas$GENE_SYMBOL






source(paste0(script_dir, 'ppmi/network_utils.R'))
############# MY ADDITION 
# 1. DE MIRS 
# 2. DE GENES 
mirnas_sig_visit

mirnas_sig_factor<-mirnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
rnas_sig_factor<-rnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)


## compare clust 1,2
mirnas_sig_factor_V08_cl2<-mirnas_sig_V08_cl2 %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
mirnas_sig_factor_V08_cl1<-mirnas_sig_V08_cl1 %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
mirnas_sig_factor_V08_cl3<-mirnas_sig_V08_cl3 %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)

rnas_sig_factor_V08_cl2<-rnas_sig_V08_cl2  %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)
rnas_sig_factor_V08_cl1<-rnas_sig_V08_cl1  %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)
rnas_sig_factor_V08_cl3<-rnas_sig_V08_cl3  %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)


intersect(rnas_sig_factor_V08_cl1$GENE_SYMBOL, rnas_sig_factor_V08_cl2$GENE_SYMBOL)


## input: 
## de_rnas: union of all
## de_mirnas: union of all clusters 
df_list = list(rnas_sig_factor_V08_cl1,rnas_sig_factor_V08_cl2,rnas_sig_factor_V08_cl3)
rnas_sig_factor_all_clusts<-Reduce(function(x, y) merge(x, y, all=TRUE), df_list)


df_list_mirs<-list(mirnas_sig_factor_V08_cl1,mirnas_sig_factor_V08_cl2,mirnas_sig_factor_V08_cl3 )
mirnas_sig_factor_all_clusts<-Reduce(function(x, y) merge(x, y, all=TRUE), df_list_mirs)

mirnas_sig_factor_all_clusts$GENE_SYMBOL
g_extended_both_clusters<-create_regulatory_net_backbone(rnas_sig_factor_all_clusts, mirnas_sig_factor_all_clusts)


#g_extended<-create_regulatory_net_backbone(rnas_sig_factor, mirnas_sig_factor)


### Now plot! 
g_extended_both_clusters
OPI_g_union<-g_extended_both_clusters
# TODO: give the layout 

set.seed(123) 
# color by logFC 
g_fc<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_BL,  de_mirnas=mirnas_BL)


g_fc_V08<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08,  de_mirnas=mirnas_V08)

g_fc_V08_cl2<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08_cl2,  de_mirnas=mirnas_V08_cl2)
g_fc_V08_cl1<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08_cl1,  de_mirnas=mirnas_V08_cl1)
g_fc_V08_cl3<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08_cl3,  de_mirnas=mirnas_V08_cl3)


g_fc_V06<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V06,  de_mirnas=mirnas_V06)
g_fc_V04<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V04,  de_mirnas=mirnas_V04)



## TODO:: save the layout 

V(g_fc)$name[is.na(V(g_fc)$color)]
V(g_fc)$group<-NA

rnas_sig_visit

## igraph plotting 

### Plot with specified coords 
#Get the coordinates of the Nodes
Coords <- layout_with_fr(g_fc_V08)# %>% 
  #as_tibble %>%
  #  bind_cols(data_frame(names = names(V(g_fc))))


Coords

sel_factor

g_fc_V06
V(g_fc_V06)$FC
g_fc_plot<-g_fc_V08_cl3
VISIT_COMP

sel_factor
V(g_fc)$group<-NA

## ADD border if significant 
# create a vertor of border colours conditional on node type
bd <- ifelse(V(g_fc_plot)$significant, "#2fc729", NA) 
bd

V(g_fc_plot)$size<-ifelse(!is.na(V(g_fc_plot)$FC), log2(abs(V(g_fc_plot)$FC)+1)*17, 0.1*20)
plot(g_fc_plot, 
vertex.color=(adjustcolor(V(g_fc_plot)$color, alpha.f = 0.8 )), 
vertex.size=V(g_fc_plot)$size,
vertex.frame.color=bd,
vertex.frame.width=3,
  layout = Coords)
title(paste('Visit: ', VISIT_COMP, ', Cluster: ', mofa_cluster_id, ', Factor: ', sel_factor, top_fr  ),cex.main=3,col.main="green")


top_fr

#V(g_fc)$shape<- ifelse(grepl("miR|hsa-let",igraph::V(g)$name), "vrectangle", "circle")
visnet <- toVisNetworkData(g_fc)
visnet$nodes$abs_FC<- abs(visnet$nodes$FC)
 

 #visualize_net(visnet)



















































