
### Load the rnas and mirnas

# 1. 

#BiocManager::install('OmnipathR')
library('OmnipathR')

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)
# interactions - proteins

mofa_cluster_id<-3
VISIT_COMP<-'V08'
load_de_by_visit_and_cluster<-function(VISIT_COMP , mofa_cluster_id, prefix='rnas_' ){
     #'
     #' @param VISIT_COMP
     #' @param mofa_cluster_id   
    rnas_visit<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c0/', VISIT_COMP ,  '/', prefix,  'de_cluster_', mofa_cluster_id, '.csv'))
    return(rnas_visit )

}

rnas_V08<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id )
rnas_sig_V08<-rnas_V08%>% dplyr::filter(mofa_sign == 'Significant')

rnas_V06<-load_de_by_visit_and_cluster(VISIT_COMP='V06', mofa_cluster_id )
rnas_V04<-load_de_by_visit_and_cluster(VISIT_COMP='V04', mofa_cluster_id )
rnas_BL<-load_de_by_visit_and_cluster(VISIT_COMP='BL', mofa_cluster_id )

rnas_sig_V06<-rnas_V06%>% dplyr::filter(mofa_sign == 'Significant')


### define the rnas fcs that will be used 
rnas_visit<-rnas_V08
rnas_sig_visit<-rnas_sig_V08



## load mirs
mirnas_V08<-load_de_by_visit_and_cluster(VISIT_COMP, mofa_cluster_id = mofa_cluster_id, prefix='mirnas_'); mirnas_V08$GENE_SYMBOL<-mirnas_V08$X
mirnas_visit<-mirnas_V08
mirnas_V06<-load_de_by_visit_and_cluster('V06', mofa_cluster_id = mofa_cluster_id, prefix='mirnas_'); mirnas_V06$GENE_SYMBOL<-mirnas_V06$X
mirnas_V04<-load_de_by_visit_and_cluster('V04', mofa_cluster_id = mofa_cluster_id, prefix='mirnas_'); mirnas_V04$GENE_SYMBOL<-mirnas_V04$X
mirnas_BL<-load_de_by_visit_and_cluster('BL', mofa_cluster_id = mofa_cluster_id, prefix='mirnas_'); mirnas_BL$GENE_SYMBOL<-mirnas_BL$X


mirnas_sig_visit<-mirnas_visit%>% dplyr::filter(mofa_sign == 'Significant')



rnas_top<-rnas_visit %>%arrange(padj) %>%  dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame
mirnas_top<-mirnas_visit %>%arrange(padj) %>% dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame

sel_factor=4; top_fr=0.05
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


mirnas_sig_factor<-mirnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
rnas_sig_factor<-rnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)










### Start loading dbs interactions, mirnas, mirtarges

## Until the DoRothEA issue gets fixed we have this here:
interactions_string <-
    import_omnipath_interactions(resources=c("SIGNOR", "STRING_talklr", "ORegAnno", "DoRothEA"))

interactions_dor <- import_transcriptional_interactions(
    resources = c("ORegAnno", "DoRothEA", "SIGNOR", "STRING_talklr" )
)


dim(interactions_dor)
dim(interactions_string)

## ----mirnatarget--------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions_mirs <-
  import_mirnatarget_interactions(resources = c("miR2Disease", "miRDeathDB", "miRTarBase", "TransmiR"))



## 1. interactions_mirs- obtain its genes 
## 2. filter by genes that have a de target OR are DE themselves 
interactions_de_mirs<-interactions_mirs %>% 
                dplyr::filter(source_genesymbol %in% mirnas_sig_factor$GENE_SYMBOL)

dim(interactions_de_mirs)
mirtargets<-interactions_de_mirs$target_genesymbol

# find out which mir targets have a de interaction
# TODO: this is only one way ie. if the mir is a source. 
interactions_of_mirtargets_tar_is_de<-interactions_string %>% dplyr::filter(
    source_genesymbol %in% mirtargets) %>% dplyr::filter(
    target_genesymbol %in% rnas_sig_factor$GENE_SYMBOL

)
dim(interactions_of_mirtargets_tar_is_de)


interactions_of_mirtargets_tar_is_de
# 2. Now do the filter of gene targets 
# 1. either de genes (in the network ) or de themselves 

interactions_de_mirs_targets_filt<-interactions_de_mirs %>%
        dplyr::filter( target_genesymbol %in% interactions_of_mirtargets_tar_is_de$source_genesymbol | 
        target_genesymbol %in% rnas_sig_factor$GENE_SYMBOL # CAREFUL to access these with gene symbol!!
         )



interactions_mirs$source_genesymbol
interactions_de_mirs_targets_filt$source_genesymbol[grep( 'miR-7-5p',interactions_de_mirs_targets_filt$source_genesymbol)]
interactions_de_mirs_targets_filt$target_genesymbol[grep( 'SNCA',interactions_de_mirs_targets_filt$target_genesymbol)]


interactions_de_mirs_targets_filt







g_mirs_targets_filt<-interaction_graph(interactions_de_mirs_targets_filt)
g_targets<-interaction_graph(interactions_of_mirtargets_tar_is_de)

g_mirs_targets_filt




############# Add also gene-gene interactions that are not in the graph already ###

interactions_dor_target_genes<-interactions_dor %>%
    dplyr::filter(target_genesymbol %in% c(rnas_sig_factor$GENE_SYMBOL ) ) %>%
     dplyr::filter( source_genesymbol %in% c(rnas_sig_factor$GENE_SYMBOL )) 


dim(interactions_dor_target_genes)
g_genes_inter<-interaction_graph(interactions_dor_target_genes)


g_extended<-union(g_mirs_targets_filt, g_targets) # mir-gene-gene interactions
#g_extended<-union(g_extended, g_genes_inter) ## add gene-gene interactions




### Now plot! 
OPI_g_union<-g_extended
# TODO: give the layout 

set.seed(123) 
# color by logFC 
g_fc<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_BL,  de_mirnas=mirnas_BL)


V(g_fc)$name[is.na(V(g_fc)$color)]
V(g_fc)$group<-NA

#V(g_fc)$shape<- ifelse(grepl("miR|hsa-let",igraph::V(g)$name), "vrectangle", "circle")
visnet <- toVisNetworkData(g_fc)
visnet$nodes$abs_FC<- abs(visnet$nodes$FC)
 

 visualize_net(visnet)





### Plot with specified coords 
#Get the coordinates of the Nodes
Coords <- layout_with_fr(g_fc) %>% 
  as_tibble %>%
    bind_cols(data_frame(names = names(V(g_fc))))

Coords
#get the coordinates of the remaining Nodes
  NetCoords <- data_frame(names = names(V(g_fc))) %>%
    left_join(Coords, by= "names")










































