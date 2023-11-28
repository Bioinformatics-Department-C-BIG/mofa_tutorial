

BiocManager::install('OmnipathR')
library('OmnipathR')


interactions <- import_omnipath_interactions( resources = c('SignaLink3', 'PhosphoSite', 'SIGNOR') )

# enzyme -PTM relationships
 enzsub <- import_omnipath_enzsub(resources = c('PhosphoSite', 'SIGNOR'))

# convert to igraph objects
ptms_g = ptms_graph(ptms = enzsub) 
OPI_g = interaction_graph(interactions = interactions)
