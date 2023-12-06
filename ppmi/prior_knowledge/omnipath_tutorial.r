
## ----fig1, fig.cap="Overview of the resources featured in OmniPath. Causal resources (including activity-flow and enzyme-substrate resources) can provide direction (*) or sign and direction (+) of interactions.", echo=FALSE----
library(knitr)
knitr::include_graphics("man/figures/page1_1.png")

## ----installation, eval=FALSE-------------------------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("OmnipathR")

## ----libraries, message=FALSE-------------------------------------------------------------------------------
library(OmnipathR)
library(tidyr)
library(dnet)
#BiocManager::install('dnet')
library(gprofiler2)

## ----interactions-------------------------------------------------------------------------------------------
## We check some of the different interaction databases
get_interaction_resources()

## The interactions are stored into a data frame.
interactions <-
    import_omnipath_interactions(resources=c("SignaLink3","PhosphoSite",
    "SIGNOR"))

## We visualize the first interactions in the data frame.
print_interactions(head(interactions))

## ----sp, message=TRUE---------------------------------------------------------------------------------------
## We transform the interactions data frame into a graph
OPI_g <- interaction_graph(interactions = interactions)

## Find and print shortest paths on the directed network between proteins
## of interest:
print_path_es(shortest_paths(OPI_g,from = "TYRO3",to = "STAT3",
    output = 'epath')$epath[[1]],OPI_g)

## Find and print all shortest paths between proteins of interest:
print_path_vs(all_shortest_paths(OPI_g,from = "DYRK2",
    to = "MAPKAPK2")$res,OPI_g)

## ----clustering, message=FALSE------------------------------------------------------------------------------
## We apply a clustering algorithm (Louvain) to group proteins in
## our network. We apply here Louvain which is fast but can only run
## on undirected graphs. Other clustering algorithms can deal with
## directed networks but with longer computational times,
## such as cluster_edge_betweenness. These cluster methods are directly
## available in the igraph package.
OPI_g_undirected <- as.undirected(OPI_g, mode=c("mutual"))
OPI_g_undirected <- simplify(OPI_g_undirected)
cl_results <- cluster_fast_greedy(OPI_g_undirected)
## We extract the cluster where a protein of interest is contained
cluster_id <- cl_results$membership[which(cl_results$names == "ERBB2")]
module_graph <- induced_subgraph(OPI_g_undirected,
    V(OPI_g)$name[which(cl_results$membership == cluster_id)])

## ----fig2, echo = FALSE, fig.cap="ERBB2 associated cluser. Subnetwork extracted from the interactions graph representing the cluster where we can find the gene *ERBB2* (yellow node)"----
## We print that cluster with its interactions.
par(mar=c(0.1,0.1,0.1,0.1))
plot(module_graph, vertex.label.color="black",vertex.frame.color="#ffffff",
    vertex.size= 15, edge.curved=.2,
    vertex.color = ifelse(igraph::V(module_graph)$name == "ERBB2","yellow",
    "#00CCFF"), edge.color="blue",edge.width=0.8)

## ----pathwayextra-------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions <-
    import_pathwayextra_interactions(resources=c("BioGRID","STRING"),
    organism = 10090)

## We select all the interactions in which Amfr gene is involved
interactions_Amfr <- dplyr::filter(interactions, source_genesymbol == "Amfr" |
    target_genesymbol == "Amfr")

## We print these interactions:
print_interactions(interactions_Amfr)

## ----kinaseextra--------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions <-
    import_kinaseextra_interactions(resources=c("PhosphoPoint",
    "PhosphoSite"), organism = 10116)

## We select the interactions in which Dpysl2 gene is a target
interactions_TargetDpysl2 <- dplyr::filter(interactions,
    target_genesymbol == "Dpysl2")

## We print these interactions:
print_interactions(interactions_TargetDpysl2)

## ----ligrecextra--------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions <- import_ligrecextra_interactions(resources=c("iTALK",
    "Baccin2019"), organism=9606)

## Receptors of the CDH1 ligand.
interactions_ADM2 <- dplyr::filter(interactions, source_genesymbol == "ADM2")

## We transform the interactions data frame into a graph
OPI_g <- interaction_graph(interactions = interactions_ADM2)

## We induce a network with these genes
Induced_Network <-  dNetInduce(g=OPI_g,
    nodes_query=as.character( V(OPI_g)$name), knn=0,
    remove.loops=FALSE, largest.comp=FALSE)

## ----fig3, echo = FALSE, fig.cap="Ligand-receptor interactions for the ADM2 ligand."------------------------
## We print the induced network
par(mar=c(0.1,0.1,0.1,0.1))
plot(Induced_Network, vertex.label.color="black",
    vertex.frame.color="#ffffff",vertex.size= 20, edge.curved=.2,
    vertex.color =
        ifelse(igraph::V(Induced_Network)$name %in% c("ADM2"),
        "yellow","#00CCFF"), edge.color="blue",edge.width=0.8)

## ----dorothea-----------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions <- import_dorothea_interactions(
    resources = c("DoRothEA"),
    dorothea_levels = 'A',
    organism = 9606
)

## Until the DoRothEA issue gets fixed we have this here:
interactions <- import_transcriptional_interactions(
    resources = c("ORegAnno", "DoRothEA")
)

## We select the most confident interactions for a given TF and we print
## the interactions to check the way it regulates its different targets
interactions_A_GLI1  <- dplyr::filter(
    interactions,
    source_genesymbol == "GLI1"
)


print_interactions(interactions_A_GLI1)

## ----mirnatarget--------------------------------------------------------------------------------------------
## We query and store the interactions into a dataframe
interactions <-
  import_mirnatarget_interactions(resources = c("miR2Disease", "miRDeathDB"))

## We select the interactions where a miRNA is interacting with the TF
## used in the previous code chunk and we print these interactions.
interactions_miRNA_GLI1 <-
    dplyr::filter(interactions, target_genesymbol == "GLI1")

print_interactions(interactions_miRNA_GLI1)

## We transform the previous selections to graphs (igraph objects)
OPI_g_1 <- interaction_graph(interactions = interactions_miRNA_GLI1)
OPI_g_2 <- interaction_graph(interactions = interactions_miRNA_GLI1)

## ----fig4, echo = FALSE, fig.cap="miRNA-TF-target network. Schematic network of the miRNA (red square nodes) targeting \textit{GLI1} (yellow node) and the genes regulated by this TF (blue round nodes)."----
## We print the union of both previous graphs
par(mar=c(0.1,0.1,0.1,0.1))
plot(OPI_g_1 %u% OPI_g_2, vertex.label.color="black",
    vertex.frame.color="#ffffff",vertex.size= 20, edge.curved=.25,
    vertex.color = ifelse(grepl("miR",igraph::V(OPI_g_1 %u% OPI_g_2)$name),
    "red",ifelse(igraph::V(OPI_g_1 %u% OPI_g_2)$name == "GLI1",
    "yellow","#00CCFF")), edge.color="blue",
    vertex.shape = ifelse(grepl("miR",igraph::V(OPI_g_1 %u% OPI_g_2)$name),
    "vrectangle","circle"),edge.width=0.8)

## ----small-molecules----------------------------------------------------------------------------------------
trametinib_interactions <- import_small_molecule_protein_interactions(
    sources = 'TRAMETINIB'
)
print_interactions(trametinib_interactions)

## ----PTMs---------------------------------------------------------------------------------------------------
## We check the different PTMs databases
get_enzsub_resources()

## We query and store the enzyme-PTM interactions into a dataframe.
## No filtering by databases in this case.
enzsub <- import_omnipath_enzsub()

## We can select and print the reactions between a specific kinase and
## a specific substrate
print_interactions(dplyr::filter(
    enzsub,
    enzyme_genesymbol == "MAP2K1",
    substrate_genesymbol == "MAPK3"
))

## In the previous results, we can see that enzyme-PTM relationships do not
## contain sign (activation/inhibition). We can generate this information
## based on the protein-protein OmniPath interaction dataset.
interactions <- import_omnipath_interactions()
enzsub <- get_signed_ptms(enzsub, interactions)

## We select again the same kinase and substrate. Now we have information
## about inhibition or activation when we print the enzyme-PTM relationships
print_interactions(dplyr::filter(enzsub,enzyme_genesymbol=="MAP2K1",
    substrate_genesymbol=="MAPK3"))

## We can also transform the enzyme-PTM relationships into a graph.
enzsub_g <- enzsub_graph(enzsub = enzsub)

## We download PTMs for mouse
enzsub <- import_omnipath_enzsub(
    resources = c("PhosphoSite", "SIGNOR"),
    organism = 10090
)

## ----complexes----------------------------------------------------------------------------------------------
## We check the different complexes databases
get_complex_resources()

## We query and store complexes from some sources into a dataframe.
complexes <- import_omnipath_complexes(resources=c("CORUM", "hu.MAP"))

## We check all the molecular complexes where a set of genes participate
query_genes <- c("WRN","PARP1")

## Complexes where any of the input genes participate
complexes_query_genes_any <- unique(get_complex_genes(complexes,query_genes,
    total_match=FALSE))

## We print the components of the different selected components
head(complexes_query_genes_any$components_genesymbols,6)

## Complexes where all the input genes participate jointly
complexes_query_genes_join <- unique(get_complex_genes(complexes,query_genes,
    total_match=TRUE))

## We print the components of the different selected components
complexes_query_genes_join$components_genesymbols

## ----enrichment---------------------------------------------------------------------------------------------
genes_complex <-
  unlist(strsplit(complexes_query_genes_join$components_genesymbols, "_"))

## We can perform an enrichment analyses with the genes in the complex
EnrichmentResults <- gost(genes_complex, significant = TRUE,
    user_threshold = 0.001, correction_method = c("fdr"),
    sources=c("GO:BP","GO:CC","GO:MF"))

## We show the most significant results
EnrichmentResults$result %>%
  dplyr::select(term_id, source, term_name,p_value) %>%
  dplyr::top_n(5,-p_value)

## ----complex_annotations------------------------------------------------------------------------------------
## We check the different annotation databases
get_annotation_resources()

## We can further investigate the features of the complex selected
## in the previous section.

## We first get the annotations of the complex itself:
annotations <- import_omnipath_annotations(proteins=paste0("COMPLEX:",
  complexes_query_genes_join$components_genesymbols))

head(dplyr::select(annotations,source,label,value),10)

## ----annotations_components---------------------------------------------------------------------------------
annotations <- import_omnipath_annotations(
    proteins = genes_complex,
    resources = "NetPath"
)

dplyr::select(annotations, genesymbol, value)

## ----subcell_loc--------------------------------------------------------------------------------------------
annotations <-import_omnipath_annotations(
    proteins = genes_complex,
    resources = "ComPPI"
)

## ----annot_spread-------------------------------------------------------------------------------------------
tidyr::spread(annotations, label, value) %>%
dplyr::arrange(desc(score)) %>%
dplyr::top_n(10, score)

## ----annot_wide---------------------------------------------------------------------------------------------
signaling_pathways <- import_omnipath_annotations(
    resources = 'SignaLink_pathway',
    wide = TRUE
)

## ----intercell----------------------------------------------------------------------------------------------
## We check some of the different intercell categories
get_intercell_generic_categories()

## We import the intercell data into a dataframe
intercell <- import_omnipath_intercell(scope = 'generic',
    aspect = 'locational')

## We check the intercell annotations for the individual components of
## our previous complex. We filter our data to print it in a good format
dplyr::filter(intercell,genesymbol %in% genes_complex) %>%
dplyr::distinct(genesymbol, parent, .keep_all = TRUE) %>%
dplyr::select(category, genesymbol, parent) %>%
dplyr::arrange(genesymbol)

## ----intercell_quality--------------------------------------------------------------------------------------
icn <- import_intercell_network(high_confidence = TRUE)

## ----intercell_filter---------------------------------------------------------------------------------------
icn <-
    import_intercell_network() %>%
    filter_intercell_network(
        min_curation_effort = 1,
        consensus_percentile = 33
    )

## ----close_dev----------------------------------------------------------------------------------------------
## We close graphical connections
while (!is.null(dev.list()))  dev.off()

## ----sessionInfo, echo=FALSE--------------------------------------------------------------------------------
sessionInfo()
