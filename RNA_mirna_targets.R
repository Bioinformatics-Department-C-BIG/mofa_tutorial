

# 2. List of mirs ####
#gsea_results_fname<-paste0(mir_results_file,'_mieaa_res.csv' )
library('igraph')

source(paste0(script_dir, '/ppmi/emap_utils.R'))

pvalueCutoff=1

mirs_de<-read.csv(paste0(outdir, '/trajectories/most_sig_over_time_2feats_TRUE_cl_fs_2-12_miRNA_kmeans_groupingde_group_1.csv'))


#mirs_de<-read.csv(paste0(outdir, '/trajectories/most_sig_vs_control_2feats_TRUE_cl_fs_2-12_miRNA_kmeans_groupingde_group_1.csv'))

mirs=gsub( '\\.','-', de_group_vs_control_and_time2)

mirs=gsub( '\\.','-', most_sig_over_time1$symbol); clust_index=1
mirs=gsub( '\\.','-', mirs_de$symbol); clust_index=1

mirs


#mirs=gsub( '\\.','-', selected_mirs)
## NOTE : the more mirs you add the less gene targets you get back as significant....

mieaa_results_fname = paste0(outdir, '/trajectories/enrichment/', clust_name, clust_index)

dir.create(mieaa_results_fname)
test_type='ORA'
if (file.exists(mieaa_results_fname)){
  ### Load enrichment results if available
  mieaa_all_gsea<-read.csv(mieaa_results_fname, header=TRUE)
  ### TODO: Rerun with updated pvalue cutoff 
}else{
  ## otherwise run GSEA analysis 
  
  
  mieaa_all_gsea <- rba_mieaa_enrich(test_set = mirs,
                                     mirna_type = "mature",
                                     test_type = test_type,
                                     species = 'Homo sapiens',
                                     sig_level=pvalueCutoff
  )
  
  
  write.csv(mieaa_all_gsea, mieaa_results_fname, row.names = FALSE)
  
}


#### 
colnames(mieaa_all_gsea)<-make.names(  colnames(mieaa_all_gsea))

mieaa_all_gsea_sig<-mieaa_all_gsea %>%
  dplyr::filter(Category %in%c('GO Biological process (miRPathDB)')) %>%
  dplyr::filter(P.adjusted<0.05)

mieaa_all_gsea_sig$Subcategory
#View(mieaa_all_gsea_sig)

mieaa_targets<-mieaa_all_gsea %>%
  dplyr::filter(Category %in%c('Target genes (miRTarBase)')) %>%
  dplyr::filter(P.value<0.05)


#write.csv(mieaa_targets, paste0(mieaa_results_fname, 'mir_targets.csv') )
#head(mieaa_targets)
mieaa_targets
to_load<-paste0(outdir,'/trajectories/enrichment/NP2PTOT_clust1mir_targets.csv')





################# 
BiocManager::install('miRTarBase')




mieaa_targets_load<-read.csv(to_load)
mieaa_targets_load

############# FILTER BY THE RELEVANT ONES ####

dim(mieaa_targets)
#View(mieaa_targets)




rnas_de<-read.csv(paste0(outdir, '/trajectories/most_sig_over_time_2feats_TRUE_cl_fs_2-12_RNA_kmeans_groupingde_group_1.csv'))
rnas_de<-read.csv(paste0(outdir, '/trajectories/most_sig_vs_control_2feats_TRUE_cl_fs_2-12_RNA_kmeans_groupingde_group_1.csv'))

ens<-as.character(rnas_de$symbol)
symb<-get_symbols_vector(ens)
rnas_de$GENE_SYMBOL<-symb




list_of_rnas = rnas_de$GENE_SYMBOL 
#list_of_rnas = de_group_vs_control1[de_group_vs_control1$symbol  %in% top_factor_feats,
#                                   'GENE_SYMBOL']

list_of_rnas
#####
selected_rnas<-intersect(list_of_rnas, mieaa_targets_load$Subcategory)
selected_rnas
mieaa_targets_de_rnas<-mieaa_targets_load %>% dplyr::filter(Subcategory %in%selected_rnas)
selected_mirs_from_targets<-mieaa_targets_de_rnas

mirna_targets_edgel<-get_long_mir_targets(mieaa_targets_de_rnas)
  
selected_mirs<-mirna_targets_edgel$variable

mirna_targets_el<-mirna_targets_edgel

g<-graph_from_edgelist(as.matrix(mirna_targets_el[,c(1,2)]), directed = TRUE)
g$layout <- layout_with_kk

p<-plot.igraph(g )
bc <- betweenness(g,v = V(g), directed = FALSE )
bc[order(bc)]

deg<-degree(g)
deg[deg>5]

bc_remove_mirs<-deg[deg<5]
bc_remove_mirs<-bc_remove_mirs[startsWith(names(bc_remove_mirs), 'hsa')]
#genes_to_keep<-V(g)$name[!startsWith(V(g)$name, 'hsa')]
deg
g_filt<-delete_vertices(g, V(g)$name %in% c(names(bc_remove_mirs)))
g_filt<-remove_subcomponents(g_filt, subcomp_min_edge = 1)[[1]]
bc <- betweenness(g_filt,v = V(g_filt), directed = FALSE )

#nodes_to_plot<-V(g_filt)
nodes_to_plot<-names(bc[bc>80])
nodes_to_plot<-c('FOXO3','SREBF2', 'ZFHX3', 'CC2D2A') 





library(visNetwork)
library(htmlwidgets)

visnet<-toVisNetworkData(g_filt)
length(visnet$x$nodes)
visnet$nodes$font.size=25
#visnet$nodes$font<-list(size=25,mod='bold' )
mst<-g_filt
mst.communities <- edge.betweenness.community(mst, weights=NULL, directed=FALSE)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <-  mst.communities$membership + 1

visnet$nodes$clusters<-mst.communities$membership

visnet$nodes$font.size=60
visnet$ed

visnet$nodes$size=log2(bc+1)*3

visnet$nodes$font.color='black'
visnet$edges$color='darkgray'
visnet$edges$width=3

visnet$nodes$type<-visnet$nodes$label
visnet$nodes$type<-ifelse(startsWith(visnet$nodes$label, 'hsa'), 'mir', 'rna')
visnet$nodes$color<-ifelse(startsWith(visnet$nodes$label, 'hsa'),'orange', 'lightblue')
#icon = list( face ='FontAwesome', code = "f0c0")


vis_net_vis<-visNetwork(visnet$nodes, visnet$edges, 
           color='color') %>%
  addFontAwesome()%>%
 # visIgraphLayout(layout = 'layout_nicely', smooth = TRUE)
 visIgraphLayout(layout = 'layout.fruchterman.reingold', smooth = TRUE)#%>%
  #visLayout(randomSeed = 150 )
#improvedLayout =TRUE,

                  #physics=TRUE)# same as   visLayout(hierarchical = TRUE) 
vis_net_vis


visSave(vis_net_vis, file = paste0(mieaa_results_fname, '.html'))

saveWidget(vis_net_vis, file = paste0(mieaa_results_fname, '.html'))





