

library('igraph')
library('visNetwork')
library(purrr)
library(dplyr)

### NEEDED to map paths to go ids 
#
lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 

### First match withgo_ids 
go_ids<-read.csv(paste0(data_dir,'ppmi/ppmi_data/go_pathway_info.txt'), sep='\t')
go_ids$GOID
go_ids$TERM[grep( 'nervous system development', go_ids$TERM )]
go_ids$TERM_standardized=standardize_go_names(go_ids$TERM)

#f_pvals_sig<-lapply(f_pvals, function(x){x[x$fish<0.05,]})







########
jaccard_from_numbers <- function(intersection, length_a, length_b) {
  union = length_a + length_b - intersection
  return (intersection/union)
}



input_gene_overlap_f<-paste0(outdir, '/enrichment/', 'GMinadakis_gene_overlap_mofa.csv' )
input_gene_overlap_f<-paste0(outdir, '/enrichment/', 'scored_nets_/NET_union_pathways0.05_V08_p_anova_FALSEpval_TRUE_GOids.csv' )
input_gene_overlap_f
#lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 




input_gene_overlap<-read.csv(input_gene_overlap_f, row.names = 1)
head(input_gene_overlap)
input_gene_overlap$jaccard<- jaccard_from_numbers(input_gene_overlap$common,input_gene_overlap$N1, input_gene_overlap$N2)
hist(log2(input_gene_overlap$jaccard))

grep( 'GO:0007399',input_gene_overlap$node2)
which(input_gene_overlap$common<1)
mofa_net<-input_gene_overlap[input_gene_overlap$common>1,]
hist(mofa_net$edgeScore);
mofa_net<-mofa_net[mofa_net$edgeScore>0.03,]

hist(log2(mofa_net$edgeScore*100))
mofa_net$node1 

emap_mofa<-graph_from_data_frame(mofa_net, directed = FALSE)


#### Set all attributes of the emap plot ####





terms<-go_ids$TERM[match(V(emap_mofa)$name,go_ids$GOID )]
node_factors<-f_pvals_merged_sig$Least_factor[match(V(emap_mofa)$name  , f_pvals_merged_sig$GOID)]
## factor info: select only from the dataframe with specific names 
### Here we merge the network with mofa attribute data
# first merge with 
all_factor_info<-f_pvals_merged_fs_sig[match(V(emap_mofa)$name  , f_pvals_merged_sig$GOID),]
all_factor_info<-cbind(all_factor_info,f_pvals_merged_sig[match(V(emap_mofa)$name  , f_pvals_merged_sig$GOID),])


node_fs<-data.frame(apply(all_factor_info[,names(sel_factors)], 2, function(x){x<0.05} ))
colnames(node_fs)<-names(sel_factors)
dim(node_fs)

V(emap_mofa)$label <- terms
V(emap_mofa)$label_description <- terms
E(emap_mofa)$my_edge_score <- E(emap_mofa)$jaccard
E(emap_mofa)$my_edge_score <- log2(E(emap_mofa)$edgeScore*100)
E(emap_mofa)$weight <- E(emap_mofa)$my_edge_score
V(emap_mofa)$groups <- node_factors


#emap_mofa <- 

colnames(f_pvals_merged_sig_fs)
for (i in 1:length(node_fs)){
  emap_mofa<-set_vertex_attr(emap_mofa, name = (colnames(node_fs)[i]), 
                               index = V(emap_mofa), value = node_fs[,i])
  
  pval_name<-paste0('pval_',(colnames(all_factor_info)[i]))
  emap_mofa<-set_vertex_attr(emap_mofa, name = pval_name, 
                             index = V(emap_mofa), value = all_factor_info[,i])
  
  rank_name<-paste0((colnames(all_factor_info)[i]), '_rank')
  
  emap_mofa<-set_vertex_attr(emap_mofa, name = rank_name, 
                             index = V(emap_mofa), value = all_factor_info[,rank_name])
  
  
  
}

V(emap_mofa)$Factor5_rank


#### Remove small components ####

#emap_mofa_rem<-remove_subcomponents(emap_mofa)


##############################
#### Convert to visnetwork #### 


choose_f<-NULL
choose_f<-names(sel_factors[1])
choose_f
### If we filter by a factor
if (is.null(choose_f)){
  emap_mofa_filt=emap_mofa
  factor_filters<-names(sel_factors)
}else{
  
  ### Filter by Factor 
  to_remove_vs<-V(emap_mofa)[!(vertex_attr(emap_mofa, choose_f)) | is.na(vertex_attr(emap_mofa, choose_f)) ]
  to_remove_vs
  #### 1. First select Factor 
  emap_mofa_filt<-delete.vertices(emap_mofa, to_remove_vs )
  factor_filters<-choose_f
}
length(emap_mofa)
length(emap_mofa_filt)
#table(components(emap_mofa_rem)$membership)



#### 2. Remove low ranks ####
high_quant_t<-0.95
high_quant_t<-1
for (factor in c(factor_filters) ){
  print(factor)
  up_quant<-quantile(vertex_attr(emap_mofa_filt,  paste0(factor, '_rank' ) ) ,high_quant_t)
  vert_gt_50<-V(emap_mofa_filt)$name[vertex_attr(emap_mofa_filt,  paste0(factor, '_rank' ) ) >up_quant]
  vert_gt_50<-vert_gt_50[!is.na(vert_gt_50)]
  vert_gt_50
  emap_mofa_filt<-delete.vertices(emap_mofa_filt,vert_gt_50)
}
length(emap_mofa_filt)
  
#### 3. Remove subcomponents of the result ####

emap_mofa_filt<-remove_subcomponents(emap_mofa_filt)
emap_mofa_filt
paste0(choose_f, '_rank' ) 


vis_emap_mofa


##########################################
#### Now convert it to a visnetwork ####

vis_emap_mofa<-toVisNetworkData(emap_mofa_filt)

### make the labeling less dense 
all_n<-length(vis_emap_mofa$nodes$label_description)
vis_emap_mofa$nodes$label_description[c(TRUE,FALSE)]<-' '


### Add labels
vis_emap_mofa$nodes$label<-vis_emap_mofa$nodes$label_description
## Set width to the edge weight
vis_emap_mofa$edges$width<-vis_emap_mofa$edges$weight
# Set length of edge to be inverse the weight
vis_emap_mofa$edges$length<-1/vis_emap_mofa$edges$weight*130
#vis_emap_mofa$edges$length<-log2((1-E(emap_mofa)$jaccard)*100)
hist(vis_emap_mofa$edges$length)

  #set.vertex.attribute(emap_mofa,



### SET NODE COLOR 
vis_emap_mofa$nodes$factor<-f_pvals_merged_sig$min_rank[match( vis_emap_mofa$nodes$id, f_pvals_merged_sig$GOID)]
which(is.na(vis_emap_mofa$nodes$factor))
unique(vis_emap_mofa$nodes$factor)


# TODO CREATE FUNCTIONS
cols_pal<-RColorBrewer::brewer.pal(4, name='Spectral') # red, 2: orange 3: green , blue
cols_pal 
grey_col<-'#808080'
## fact:1: red 3: orange , 4: green, 14: blue 
# pvalue<- factor  1 , 3,  4, 14 

mapdf <- data.frame(old=c(names(sel_factors)),new=cols_pal[1:length(sel_factors)])
vis_emap_mofa$nodes$groups
### COLOR BY FACTOR
if (is.null(choose_f)){
  vis_emap_mofa$nodes$color<-mapdf$new[match(vis_emap_mofa$nodes$factor, mapdf$old)]
  #vis_emap_mofa$nodes$color<-vis_emap_mofa$nodes[, choose_f]
  #vis_emap_mofa$nodes$color<-ifelse(vis_emap_mofa$nodes[, choose_f], cols_pal[1], cols_pal[2])
  #vis_emap_mofa$nodes$color[is.na(vis_emap_mofa$nodes[, choose_f])]<-grey_col
}else{
  ### COLOR BY P-VALUE, FILTER BY P-VALUE
  # Create continuous color palette.
  eigScalePal <- colorRampPalette(c('#E0F4FF','#003049'))
  
  # Match palette to centrality vector.
  un_ps<-length(vis_emap_mofa$nodes[,  paste0('pval_',choose_f )])
  vis_emap_mofa$nodes$color <- eigScalePal(un_ps)[cut(vis_emap_mofa$nodes[,  paste0('pval_',choose_f )], breaks=un_ps)]
  
}
choose_f





unique(vis_emap_mofa$nodes$groups)
unique(vis_emap_mofa$nodes$label)
vis_emap_mofa$nodes$label




choose_f

visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges, 
           height = "1000px", width='2000px', main=paste0(choose_f), 
          submain='Color by pvalue')%>%
  visEdges(width=width)%>%
 #visOptions(selectedBy= list(variable="label",multiple=T)) %>%
  visNodes(font=list(size=30) )%>%
  visIgraphLayout(layout = 'layout.fruchterman.reingold') %>%# same as   visLayout(hierarchical = TRUE) 
  
  visSave(file = paste0(outdir, '/enrichment/all_factors',pval_to_use, choose_f,'.html'))

vis_emap_mofa$nodes$groups[is.na(vis_emap_mofa$nodes$groups)]<-'not_in_mofa'


visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges, 
           main=paste0(choose_f))%>%
  #visOptions(selectedBy= list(variable="group",multiple=T)) %>%
  visEdges(width=width, 
           length=length)%>%
  visIgraphLayout(layout = 'layout_with_drf') %>%# same as   visLayout(hierarchical = TRUE) 
visNodes(font=list(size=50) )
#%>%
#  visGroups(groupname = "fish.x", color = "red", shape = "triangle")%>%
#  visGroups(groupname = "fish.y.y", color = "yellow", shape = "triangle")             


                       
                       