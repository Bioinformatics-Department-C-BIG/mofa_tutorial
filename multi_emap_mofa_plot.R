

library('igraph')
library('visNetwork')
### NEEDED to map paths to go ids 

lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 

go_ids<-read.csv(paste0(data_dir,'ppmi/ppmi_data/go_pathway_info.txt'), sep='\t')
go_ids$GOID
go_ids$TERM[grep( 'nervous system development', go_ids$TERM )]
go_ids$TERM_standardized=standardize_go_names(go_ids$TERM)

#f_pvals_sig<-lapply(f_pvals, function(x){x[x$fish<0.05,]})
f_pvals_merged<-f_pvals %>% reduce(full_join , by='Description')

dim(f_pvals_merged)

f_pvals_merged 

f_pvals_merged_fs<-f_pvals_merged[,c('fish.x', 'fish.x.x', 'fish.y', 'fish.y.y')]
rownames(f_pvals_merged_fs)<-f_pvals_merged$Description

### Apply Ranking: 

ranks_fs<-as.data.frame(apply(f_pvals_merged_fs, 2, rank))
tot_rank=rowSums(ranks_fs)
min_rank<-colnames(ranks_fs)[apply(ranks_fs,1,which.min)]
filter_paths_gt_top_n<-ranks_fs

f_pvals_merged_fs$Least_value<-rowMin(as.matrix(f_pvals_merged_fs))
f_pvals_merged_fs$Least_factor<-colnames(f_pvals_merged_fs)[apply(f_pvals_merged_fs,1,which.min)]
f_pvals_merged_fs$tot_rank=tot_rank
f_pvals_merged_fs$min_rank=min_rank

# filter
f_pvals_merged_fs<-f_pvals_merged_fs[f_pvals_merged_fs$Least_value<0.05,]
dim(f_pvals_merged_fs)
###

f_pvals_merged_fs$Least_factor
rownames(f_pvals_merged_fs)[grep( 'nervous system',rownames(f_pvals_merged_fs) )]
# rename descriptions 

rownames(f_pvals_merged_fs)<-standardize_go_names(rownames(f_pvals_merged_fs))
f_pvals_merged_fs$GOID<-go_ids$GOID[match(rownames(f_pvals_merged_fs),go_ids$TERM_standardized )]
f_pvals_merged_fs[, c('GOID','Least_factor')]


table(f_pvals_merged_fs$Least_factor)







########
jaccard_from_numbers <- function(intersection, length_a, length_b) {
  union = length_a + length_b - intersection
  return (intersection/union)
}

input_gene_overlap_f<-paste0(outdir, '/enrichment/', 'GMinadakis_gene_overlap_mofa.csv' )

#lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 




input_gene_overlap<-read.csv(input_gene_overlap_f, row.names = 1)
head(input_gene_overlap)
input_gene_overlap$jaccard<- jaccard_from_numbers(input_gene_overlap$common,input_gene_overlap$N1, input_gene_overlap$N2)
hist(log2(input_gene_overlap$jaccard))

grep( 'GO:0007399',input_gene_overlap$node2)
which(input_gene_overlap$common<1)
mofa_net<-input_gene_overlap[input_gene_overlap$common>1,]
hist(mofa_net$edgeScore);
mofa_net<-mofa_net[mofa_net$edgeScore>0.08,]

hist(log2(mofa_net$edgeScore*100))
mofa_net$node1 

emap_mofa<-graph_from_data_frame(mofa_net, directed = FALSE)
V(emap_mofa)$name





nodes<-V(emap_mofa)



terms<-go_ids$TERM[match(V(emap_mofa)$name,go_ids$GOID )]
node_factors<-f_pvals_merged_fs$Least_factor[match(V(emap_mofa)$name  , f_pvals_merged_fs$GOID)]
node_factors

V(emap_mofa)$label <- terms

#E(emap_mofa)$my_edge_score <- E(emap_mofa)$edgeScore
E(emap_mofa)$jaccard
head(E(emap_mofa))
E(emap_mofa)$my_edge_score <- E(emap_mofa)$jaccard
E(emap_mofa)$my_edge_score <- log2(E(emap_mofa)$edgeScore*100)
E(emap_mofa)$my_edge_score

E(emap_mofa)$weight <- E(emap_mofa)$my_edge_score

V(emap_mofa)$groups <- node_factors

plot(emap_mofa, edge.width=E(emap_mofa)$weight)


vis_emap_mofa<-toVisNetworkData(emap_mofa)
merged_results$Description
vis_emap_mofa$nodes$label
#vis_emap_mofa$nodes$pvalue<-
vis_emap_mofa$nodes$label<-terms
vis_emap_mofa$edges$width<-vis_emap_mofa$edges$weight
vis_emap_mofa$edges$width
vis_emap_mofa$edges$length<-1/vis_emap_mofa$edges$weight*100

#vis_emap_mofa$edges$length<-log2((1-E(emap_mofa)$jaccard)*100)

hist(vis_emap_mofa$edges$length)

  #set.vertex.attribute(emap_mofa,



### SET NODE COLOR 
f_pvals_merged_fs[, c('GOID','Least_factor')]



vis_emap_mofa$nodes$factor<-f_pvals_merged_fs$Least_factor[match( vis_emap_mofa$nodes$id, f_pvals_merged_fs$GOID)]
vis_emap_mofa$nodes$factor<-f_pvals_merged_fs$min_rank[match( vis_emap_mofa$nodes$id, f_pvals_merged_fs$GOID)]

which(is.na(vis_emap_mofa$nodes$factor))
unique(vis_emap_mofa$nodes$factor)



cols_pal<-RColorBrewer::brewer.pal(4, name='Spectral') # red, 2: orange 3: green , blue
cols_pal 
## fact:1: red 3: orange , 4: green, 14: blue 
# pvalue<- factor  1 , 3,  4, 14 

mapdf <- data.frame(old=c("fish.x"  , "fish.y" ,"fish.x.x", 'fish.y.y'),new=cols_pal)
vis_emap_mofa$nodes$color<-mapdf$new[match(vis_emap_mofa$nodes$groups, mapdf$old)]
#vis_emap_mofa$nodes$group<-vis_emap_mofa$nodes$color
unique(vis_emap_mofa$nodes$groups)
unique(vis_emap_mofa$nodes$label)


visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges, 
           height = "1000px", width='2000px')%>%
  visEdges(width=width)%>%
 visOptions(selectedBy= list(variable="label",multiple=T)) %>%
  visNodes(font=list(size=50) )%>%

  visSave(file = paste0(outdir, '/enrichment/all_factors',pval_to_use, '.html'))


visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges)%>%
  visOptions(selectedBy= list(variable="label",multiple=T)) %>%
  visEdges(width=width, 
           length=length)%>%
  visNodes(font=list(size=50) )
#%>%
#  visGroups(groupname = "fish.x", color = "red", shape = "triangle")%>%
#  visGroups(groupname = "fish.y.y", color = "yellow", shape = "triangle")             


                       
                       