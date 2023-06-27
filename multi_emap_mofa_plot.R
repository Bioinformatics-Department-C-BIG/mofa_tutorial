

library('igraph')
library('visNetwork')
### NEEDED to map paths to go ids 
go_ids<-read.csv(paste0(data_dir,'ppmi/ppmi_data/go_pathway_info.txt'), sep='\t')
go_ids$GOID
go_ids$TERM[grep( 'nervous system development', go_ids$TERM )]
go_ids$TERM_standardized=standardize_go_names(go_ids$TERM)


f_pvals_merged<-f_pvals %>% reduce(left_join , by='Description')
dim(f_pvals_merged)

f_pvals_merged 
f_pvals_merged_fs<-f_pvals_merged[,c('fish.x', 'fish.x.x', 'fish.y', 'fish.y.y')]
rownames(f_pvals_merged_fs)<-f_pvals_merged$Description
f_pvals_merged_fs$Least_factor<-colnames(f_pvals_merged_fs)[apply(f_pvals_merged_fs,1,which.min)]
f_pvals_merged_fs$Least_factor
rownames(f_pvals_merged_fs)[grep( 'nervous system',rownames(f_pvals_merged_fs) )]
# rename descriptions 

rownames(f_pvals_merged_fs)<-standardize_go_names(rownames(f_pvals_merged_fs))
f_pvals_merged_fs$GOID<-go_ids$GOID[match(rownames(f_pvals_merged_fs),go_ids$TERM_standardized )]
f_pvals_merged_fs[, c('GOID','Least_factor')]
########

input_gene_overlap_f<-paste0(outdir, '/enrichment/', 'GMinadakis_gene_overlap_mofa.csv' )

#lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 

input_gene_overlap<-read.csv(input_gene_overlap_f, row.names = 1)
grep( 'GO:0007399',input_gene_overlap$node2)
which(input_gene_overlap$common<2)
mofa_net<-input_gene_overlap[input_gene_overlap$common>2,]
hist(mofa_net$edgeScore);
mofa_net<-mofa_net[mofa_net$edgeScore>0.05,]

hist(log2(mofa_net$edgeScore*100))
mofa_net$node1 

emap_mofa<-graph_from_data_frame(mofa_net, directed = FALSE)
emap_mofa





nodes<-V(emap_mofa)



terms<-go_ids$TERM[match(V(emap_mofa)$name,go_ids$GOID )]
node_factors<-f_pvals_merged_fs$Least_factor[match(V(emap_mofa)$name  , f_pvals_merged_fs$GOID)]
node_factors

V(emap_mofa)$label <- terms
E(emap_mofa)$weight <- log2(E(emap_mofa)$edgeScore*100)

V(emap_mofa)$groups <- node_factors

plot(emap_mofa, edge.width=E(emap_mofa)$weight)


vis_emap_mofa<-toVisNetworkData(emap_mofa)
merged_results$Description
vis_emap_mofa$nodes$label
#vis_emap_mofa$nodes$pvalue<-
vis_emap_mofa$nodes$label<-terms
vis_emap_mofa$edges$width<-vis_emap_mofa$edges$weight
vis_emap_mofa$edges$width


  #set.vertex.attribute(emap_mofa,



### SET NODE COLOR 
f_pvals_merged_fs[, c('GOID','Least_factor')]



vis_emap_mofa$nodes$factor<-f_pvals_merged_fs$Least_factor[match( vis_emap_mofa$nodes$id, f_pvals_merged_fs$GOID)]
which(is.na(vis_emap_mofa$nodes$factor))
unique(vis_emap_mofa$nodes$factor)

unique(f_pvals_merged_fs$Least_factor)
cols_pal<-RColorBrewer::brewer.pal(4, name='Spectral') # 2: orange 3: green 
## fact:1: red 3: orange , 4: green, 14: blue 
# pvalue<- factor 14, 1 , 3,  4

mapdf <- data.frame(old=c("fish.x"  , "fish.y.y" ,"fish.y", 'fish.x.x'),new=cols_pal)
vis_emap_mofa$nodes$color<-mapdf$new[match(vis_emap_mofa$nodes$groups, mapdf$old)]
#vis_emap_mofa$nodes$group<-vis_emap_mofa$nodes$color
unique(vis_emap_mofa$nodes$groups)


visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges)%>%
  visEdges(width=width)%>%
 visOptions(selectedBy= list(variable="label",multiple=T)) %>%
  visSave(file = paste0(outdir, '/enrichment/all_factors',pval_to_use, '.html'))


visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges)%>%
  visOptions(selectedBy= list(variable="label",multiple=T)) %>%

  visEdges(width=width)
 # visNodes(color=group) 
#%>%
#  visGroups(groupname = "fish.x", color = "red", shape = "triangle")%>%
#  visGroups(groupname = "fish.y.y", color = "yellow", shape = "triangle")             


                       
                       