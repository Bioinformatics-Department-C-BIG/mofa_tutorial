## ----message=FALSE, warning=FALSE---------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(OmnipathR)
library(igraph)
library(ggraph)
library(magrittr)

#BiocManager::install('OmnipathR')
library('OmnipathR')

# interactions <- import_omnipath_interactions( resources = c('SignaLink3', 'PhosphoSite', 'SIGNOR') )
#interactions <- import_omnipath_interactions( resources = c('SignaLink3', 'PhosphoSite', 'SIGNOR','miRTarBase') )
interactions = import_omnipath_interactions(
  resources = c('SignaLink3'),
  organism = 9606
)

# enzyme -PTM relationships
 enzsub <- import_omnipath_enzsub(resources = c('PhosphoSite', 'SIGNOR'))

# convert to igraph objects
ptms_g = ptms_graph(ptms = enzsub) 
OPI_g = interaction_graph(interactions = interactions)


BiocManager::install("dorothea")
BiocManager::install("decoupleR")

### mirnas-rnas targets 
import_mirnatarget_interactions(
  resources = NULL,
  organism = 9606,
  fields = NULL,
  default_fields = TRUE,
  references_by_resource = TRUE,
  exclude = NULL
)





