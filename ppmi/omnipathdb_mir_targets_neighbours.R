
### Load the rnas and mirnas


#BiocManager::install('OmnipathR')
library('OmnipathR')

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)
# interactions - proteins

mofa_cluster_id<-1
VISIT_COMP<-'V08'
rnas_V08<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c0/','V08' ,  '/rnas_de_cluster_', mofa_cluster_id, '.csv'))
rnas_sig_V08<-rnas_V08%>% dplyr::filter(mofa_sign == 'Significant')

vis2='V08'
rnas_2<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c0/',vis2 ,  '/rnas_de_cluster_', mofa_cluster_id, '.csv'))
rnas_sig_2<-rnas_2%>% dplyr::filter(mofa_sign == 'Significant')

intersect(rnas_sig_2$GENE_SYMBOL, rnas_sig_V08$GENE_SYMBOL)

rnas<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c0/',VISIT_COMP ,  '/rnas_de_cluster_', mofa_cluster_id, '.csv'))


dim(rnas)

rnas$mofa_sign
rnas$padj
rnas_sig<-rnas%>% dplyr::filter(mofa_sign == 'Significant')
dim(rnas_sig)
rnas_sig$padj



mirnas<-read.csv(paste0(outdir,'/clustering/NP2PTOT_LOG_clust/3/TRUE/de_c0/', VISIT_COMP ,  '/mirnas_de_cluster_',mofa_cluster_id, '.csv'))
mirnas
mirnas$GENE_SYMBOL<-mirnas$X
mirnas
mirnas_sig<-mirnas%>% dplyr::filter(mofa_sign == 'Significant')
mirnas_sig

mirnas$padj

rnas_top<-rnas %>%arrange(padj) %>%  dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame
mirnas_top<-mirnas %>%arrange(padj) %>% dplyr::filter(mofa_sign == 'Significant') %>% as.data.frame

sel_factor=8; top_fr=0.2
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


mirnas_sig_factor<-mirnas_sig %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
rnas_sig_factor<-rnas_sig %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)










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

V(g_extended)



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


# color by logFC 
g_fc<-get_logFC_by_node(OPI_g_union)


V(g_fc)$name[is.na(V(g_fc)$color)]
V(g_fc)$group<-NA

#V(g_fc)$shape<- ifelse(grepl("miR|hsa-let",igraph::V(g)$name), "vrectangle", "circle")
visnet <- toVisNetworkData(g_fc)
visnet$nodes$abs_FC<- abs(visnet$nodes$FC)
 

 visualize_net(visnet)







































