

#BiocManager::install('OmnipathR')
library('OmnipathR')

library(OmnipathR)
library(tidyr)
library(dnet)
library(gprofiler2)
# interactions - proteins
outdir<-'/Volumes/GoogleDrive/Other computers/My computer (1) (1)/ppmi/plots/p_V08_CSF_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.2_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_V08_TRUEruv_1_split_FALSE (1)/'
mofa_cluster_id<-2

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

top_mirnas_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'miRNA',  factors=sel_factor, top_fr = top_fr))

length(top_genes_factor8)
## 1. Colour if in TOP 1% of factor
## 2. Colour if in specific pathways 
## 
top_genes_factor8

dim(rnas_top); dim(mirnas_top)
top_g<-1000
top_m<-200

top_rnas<-rnas_top[1:top_g,]
top_mirnas<-mirnas_top[1:top_m,]

top_mirnas$GENE_SYMBOL<-top_mirnas$X
top_mirnas$GENE_SYMBOL
library(dplyr)
mutate
interactions <- import_omnipath_interactions( resources = c( 'SIGNOR','STRING_talklr' ) )


top_mirnas
## FILTER THE REFERENCE





## ----dorothea-----------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions <- import_dorothea_interactions(
    resources = c("DoRothEA"),
    dorothea_levels = 'A',
    organism = 9606
)

## Until the DoRothEA issue gets fixed we have this here:
interactions_string <-
    import_omnipath_interactions(resources=c("SIGNOR", "STRING_talklr"))


interactions_dor <- import_transcriptional_interactions(
    resources = c("ORegAnno", "DoRothEA", "SIGNOR")
)


## ----mirnatarget--------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions_mirs <-
  import_mirnatarget_interactions(resources = c("miR2Disease", "miRDeathDB", "miRTarBase", "TransmiR"))


source(paste0(script_dir, 'ppmi/network_utils.R'))
############# MY ADDITION 
# 1. DE MIRS 
# 2. DE GENES 

mirnas_sig$GENE_SYMBOL
interactions_mirs
mirnas_sig_factor<-mirnas_sig %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
rnas_sig_factor<-rnas_sig %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)


# TODO: add top mirs in factor 8? 
## We select the interactions where a de miRNA is interacting with a de gene? 
interactions_de_mirs_de_genes <-interactions_mirs %>% 
   dplyr::filter( source_genesymbol %in% c(mirnas_sig_factor$GENE_SYMBOL)) %>% # mirs should be de
    dplyr::filter(target_genesymbol %in% rnas_sig_factor$GENE_SYMBOL)



# if a mir target is a tf bring its target genes  - OR
interactions_string_mirtars<-interactions_string %>%
 #   dplyr::filter( target_genesymbol %in% c(interactions_de_mirs_de_genes$target_genesymbol)) 
    dplyr::filter( source_genesymbol %in% c(interactions_de_mirs_de_genes$target_genesymbol)) 

print_interactions(interactions_string_mirtars)


dim(interactions_de_mirs_de_genes)

interactions_dor_target_genes<-interactions_dor %>%
    dplyr::filter(target_genesymbol %in% c(rnas_sig_factor$GENE_SYMBOL ) ) %>%
     dplyr::filter( source_genesymbol %in% c(rnas_sig_factor$GENE_SYMBOL )) 

dim(interactions_dor_target_genes)


interactions_de_mirs_de_genes$target_genesymbol
#interactions_dor_target_genes<-interactions_dor %>%
#    dplyr::filter( target_genesymbol %in%  interactions_de_mirs_de_genes$target_genesymbol) %>%
#        dplyr::filter( source_genesymbol %in% interactions_de_mirs_de_genes$target_genesymbol) 



interactions_de_mirs_de_genes$target_genesymbol

## TODO: Function filter or colour network by 
# 1. factor 
# 2. pathways 
# 2 data driven anticorelation
length(interactions_de_mirs_de_genes)
OPI_g_de_mirs_de_genes <- interaction_graph(interactions = interactions_de_mirs_de_genes)
OPI_g_dor_target_genes <-interaction_graph(interactions= interactions_dor_target_genes )
# TODO: add tfs with a different symbol ? 
# add gene interaction with different colour 
# add mir-gene target with a different colour
# 


OPI_g_de_mirs_de_genes
OPI_g_dor_target_genes

OPI_g_de_mirs_de_genes_targets<- union(OPI_g_dor_target_genes, OPI_g_de_mirs_de_genes)
OPI_g_de_mirs_de_genes_targets


g_fc<-get_logFC_by_node(OPI_g_de_mirs_de_genes_targets)
V(g_fc)$name
V(g_fc)$name[is.na(V(g_fc)$color)]
V(g_fc)$group<-NA

V(g_fc)$name[!(V(g_fc)$name %in% rnas_sig$GENE_SYMBOL)]
V(g_fc)$name[!(V(g_fc)$name %in% mirnas_sig$GENE_SYMBOL)]

 visnet$nodes
#V(g_fc)$shape<- ifelse(grepl("miR|hsa-let",igraph::V(g)$name), "vrectangle", "circle")
visnet <- toVisNetworkData(g_fc)


  # visnet rectangle 
   visnet$nodes$font.size=35
   visnet$nodes$size=10

    names(visnet$edges)
    visnet$edges<-visnet$edges[c('from', 'to')]
    vis_net_vis<-visNetwork(visnet$nodes, visnet$edges) %>%
               # visNodes( color =visnet$nodes$color  ) %>%
                visEdges(color='gray')

    vis_net_vis

    dir.create(paste0(outdir, '/networks/'))
    net_name=paste0('mirs_genes_', mofa_cluster_id, '_f',sel_factor,top_fr )
    net_name
    visSave(vis_net_vis, file = paste0(outdir, '/networks/',  net_name, '.html'))



## We select the most confident interactions for a given TF and we print
## the interactions to check the way it regulates its different targets
interactions_dor_filt  <- dplyr::filter(
    interactions_dor,
    source_genesymbol %in% rnas_sig$GENE_SYMBOL
)
rnas_sig
interactions_dor_filt



## ----mirnatarget--------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions_mirs <-
  import_mirnatarget_interactions(resources = c("miR2Disease", "miRDeathDB", "miRTarBase"))


## We select the interactions where a miRNA is interacting with the TF
## used in the previous code chunk and we print these interactions.
interactions_mirs_filt <-
    dplyr::filter(interactions_mirs, target_genesymbol %in% interactions_dor_filt$source_genesymbol)




print_interactions(interactions_mirs_filt)
OPI_g_2<-interaction_graph(interactions_mirs_filt)

OPI_g_2
interactions_mirs_filt$target_genesymbol
interactions_filt_genes$source_genesymbol

# select the interactors of the mir targets
mir_targets_interactions<-interactions_filt_genes %>% 
                    dplyr::filter(source_genesymbol %in% interactions_mirs_filt$target_genesymbol  |
                                  target_genesymbol %in%interactions_mirs_filt$target_genesymbol )

OPI_g_mirtargets_interactions<-interaction_graph(mir_targets_interactions)



which(rnas_sig$GENE_SYMBOL %in% c('IL4'))
 
which(top100_rnas$GENE_SYMBOL %in% mir_targets_interactions$target_genesymbol)
## ----fig4, echo = FALSE, fig.cap="miRNA-TF-target network. Schematic network of the miRNA (red square nodes) targeting \textit{GLI1} (yellow node) and the genes regulated by this TF (blue round nodes)."----
## We print the union of both previous graphs
par(mar=c(0.1,0.1,0.1,0.1))
OPI_g_1=OPI_g_mirtargets_interactions
top100_rnas$GENE_SYMBOL

top100_rnas$log2FoldChange

all_nodes<-igraph::V(OPI_g_1 %u% OPI_g_2)$name
all_nodes[grepl("miR",igraph::V(OPI_g_1 %u% OPI_g_2)$name)]


g<-union(OPI_g_1 %u% OPI_g_2)
g_names<-V(OPI_g_1 %u% OPI_g_2)$name

## fetch the logFC 
g<-get_log







V(g)$color[is.na(V(g)$group) ]<-"gray"
V(g)$color

V(g)$color = ifelse( igraph::V(g)$group == 'up', "red",
    # specified de genes are yellow
    ifelse( igraph::V(g)$group == 'down', "yellow",
    "#00CCFF"))
V(g)$color


plot(g, vertex.label.color="black",
    vertex.frame.color="#ffffff",vertex.size= 10, edge.curved=.25,
    # mirs are red 
   

    vertex.shape = ifelse(grepl("miR",igraph::V(g)$name),
    "vrectangle","circle"),edge.width=0.8)

### PLOT VISNET 
visnet<-toVisNetworkData(g)
visnet$nodes$shape
colnames(visnet$nodes)

#visnet$nodes$shape<-ifelse(grepl("miR", visnet$nodes$id),
#    "rectangle","circle")

dir.create('networks/')
library(htmlwidgets)



sys_pandoc <- find_program("pandoc")
sources <- c(Sys.getenv("RSTUDIO_PANDOC"), if (nzchar(sys_pandoc)) dirname(sys_pandoc))

visualize_net(visnet)
library('visNetwork')





plot(OPI_g_1 %u% OPI_g_2, vertex.label.color="black",
    vertex.frame.color="#ffffff",vertex.size= 20, edge.curved=.25,
    # DE mirs are 
    vertex.color = ifelse(grepl("miR",igraph::V(OPI_g_1 %u% OPI_g_2)$name), "red",
    # specified de genes are yellow
    ifelse(igraph::V(OPI_g_1 %u% OPI_g_2)$name %in%  top100_rnas$GENE_SYMBOL, "yellow",
    

    # MIRNAS are square
    "#00CCFF")), edge.color="blue",
    vertex.shape = ifelse(grepl("miR",igraph::V(OPI_g_1 %u% OPI_g_2)$name),
    "vrectangle","circle"),edge.width=0.8)













# Filter the top mirs and genes 
interactions_post_trans_filt<-interactions_post_transcriptional %>% 
        dplyr::filter(source_genesymbol %in% top_mirnas$GENE_SYMBOL |
        target_genesymbol %in% top100_rnas$GENE_SYMBOL ) 

interactions_post_trans_filt


get_interaction_resources()

interactions <- import_omnipath_interactions( resources = c('SignaLink3', 'PhosphoSite', 'SIGNOR') )
interactions <- import_omnipath_interactions( resources = c( 'SIGNOR') )

interactions

get_interaction_resources()

# enzyme-PTM relationships
enzsub <- import_omnipath_enzsub(resources = c('PhosphoSite', 'SIGNOR'))

enzsub
# convert to igraph objects
ptms_g = ptms_graph(ptms = enzsub) 
OPI_g = interaction_graph(interactions = interactions)

# 1. mirtar base - mir-gene anticorelations - check anticorelation within the cca? 
# 2. mirnas-gene 
# 3. 

# 1. filter mir-gene
# 2. 







