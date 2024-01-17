

library(MOFAdata)
data("CLL_data")
script_dir
wd<-'/Volumes/GoogleDrive/Other computers/My computer (1) (1)/ppmi/tutorials/'
paste0(wd,'/../Comparisons/Drug-compound cll.txt')
drugs_comps<-read.csv(paste0(script_dir,'/Comparisons/Drug-compound cll.txt'), sep='\t')
drugs_comps
dir.create(paste0(wd, '/data/'))
rownames(CLL_data$mRNA)<-get_symbols_vector(rownames(CLL_data$mRNA))
write.csv(CLL_data$mRNA,paste0(wd, '/data/', 'cll_rna.csv'))

get_highly_variable_matrix()
colnames(CLL_data$Methylation)

write.csv(CLL_data$Mutations,paste0(wd, '/data/', 'cll_mut.csv'))
#se<-DESeqDataSetFromMatrix(CLL_data$mRNA, colData = CLL_metadata)

hist(log2(CLL_data$mRNA))
dim(CLL_data$Drugs)
na.omit(CLL_data$mRNA_top)
CLL_data$mRNA_log<-log2(CLL_data$mRNA)
CLL_data$mRNA_top<-selectMostVariable(CLL_data$mRNA, 0.1)

dim(CLL_data$mRNA_top)
df<-CLL_data$mRNA_top
CLL_data$mRNA_top<-df[ , colSums(is.na(df))==0]
head(CLL_data$mRNA_top)
rownames(CLL_data$mRNA_top)<-get_symbols_vector(rownames(CLL_data$mRNA_top))

write.csv(CLL_data$mRNA_top,paste0(wd, '/data/xmwas/', 'cll_rna_top.csv'))
write.csv(t(CLL_data$mRNA_top),paste0(wd, '/data/', 'cll_rna_top_t.csv'),row.names = FALSE)

CLL_metadata$Class<-as.factor(CLL_metadata$IGHV)
colnames(CLL_data$mRNA_top)
CLL_metadata_sub<-CLL_metadata[match(colnames(CLL_data$mRNA_top),CLL_metadata$sample),]
dim(CLL_metadata_sub)
dim(CLL_data$mRNA_top)
dim(CLL_data$Methylation)

write.csv(CLL_metadata_sub,paste0(wd, '/data/', 'cll_metadata.csv'))
write.csv(CLL_metadata_sub,paste0(wd, '/data/xmwas/', 'cll_metadata.csv'))


CLL_data$Methylation_sub<-CLL_data$Methylation[,match(colnames(CLL_data$mRNA_top),colnames(CLL_data$Methylation))]
CLL_data$Methylation_sub_top<-selectMostVariable(CLL_data$Methylation_sub, 0.03)
dim(CLL_data$Methylation_sub_top)

write.csv(CLL_data$Methylation_sub_top,paste0(wd, '/data/xmwas/', 'cll_meth.csv'))
write.csv(t(CLL_data$Methylation_sub_top),paste0(wd, '/data/', 'cll_meth_t.csv'),row.names = FALSE)


CLL_data$Drugs_sub<-CLL_data$Drugs[,match(colnames(CLL_data$mRNA_top),colnames(CLL_data$Drugs))]
dim(CLL_data$Drugs_sub)
write.csv(CLL_data$Drugs_sub,paste0(wd, '/data/xmwas/', 'cll_drugs.csv') )
write.csv(t(CLL_data$Drugs_sub),paste0(wd, '/data/', 'cll_drugs_t.csv'),row.names = FALSE)




### XMWAS 
out_res<-paste0(wd, 'xmwas/','xmwasresults20231204144622/')
link_mat<-read.csv(paste0(out_res, 'Multidata_Network_threshold0.4_linkmatrix.txt' ), sep='\t')
name_map<-read.csv(paste0(out_res, 'NodeID_Name_mapping.txt'), sep='\t')
name_map$Name[match( link_mat$from, name_map$Node)]

drug_prefix<-link_mat$to2[startsWith(x = link_mat$to2, prefix = 'D_')]
drug_prefix2<- gsub('_[^_]*$','',drug_prefix) # remove everything after last underscore to matcgdruyg compiund
drug_prefix2
drug_prefix2_compound<-drugs_comps[ match(drug_prefix2, drugs_comps$Drug),'Compound']

which(drugs_in_an)
drugs_in_an
drugs_comps[drugs_comps$Drug =='D_078',]
name_map$Node
link_mat$from
library(igraph)
thresh<-0.63
link_mat$from2<-name_map$Name[match( link_mat$from, name_map$Node)]
link_mat$to2<-name_map$Name[match( link_mat$to, name_map$Node)]
link_mat$to2[startsWith(x = link_mat$to2, prefix = 'D_')]<-drug_prefix2_compound
link_mat_filt<-link_mat[abs(link_mat$weight)>thresh,]
dim(link_mat_filt)



library(visNetwork)
rel_net<-igraph::graph_from_edgelist(as.matrix(link_mat_filt[,4:5]), directed = TRUE)

E(rel_net)$weight <- abs(link_mat_filt$weight)
colPal <- c("#ECECEC", "#F9D597", "#EE9C77", 
            "#46887C", "#4270A4", "#786696", 
            "#A8534C", "#624E4D", "#232323")



### create comunities ####





#### FORMATTING ####
visnet$nodes$font.size=45
mst<-rel_net

mst.communities <- edge.betweenness.community(mst)
mst.clustering <- make_clusters(mst, membership=mst.communities$membership)
V(mst)$color <- mst.communities$membership + 1


visnet<-visNetwork::toVisNetworkData(rel_net)

plot(rel_net)
visnet$nodes$type[visnet$nodes$id %in% rownames(CLL_data$Methylation)]='methylation'
visnet$nodes$type[visnet$nodes$id %in% rownames(CLL_data$mRNA_top)]='mRNA'
visnet$nodes$type[visnet$nodes$id %in% drugs_comps$Compound]='drug'
visnet$nodes$type
visnet$nodes$color[visnet$nodes$type=='drug']<-"#EE9C77"
visnet$nodes$color[visnet$nodes$type=='mRNA']<-"#46887C"
visnet$nodes$color[visnet$nodes$type=='methylation']<-"#A8534C"




### Label with the type of omics


vis_net_vis$edges$weight <- link_mat_filt$weight
vis_net_vis$edges$value <- link_mat_filt$weight
visnet$edges$width=visnet$edges$weight*2


vis_net_vis<-visNetwork(visnet$nodes, visnet$edges
                        )%>%

  addFontAwesome()%>%
  # visIgraphLayout(layout = 'layout_nicely', smooth = TRUE)
  #visIgraphLayout(layout = 'layout.fruchterman.reingold', smooth = TRUE)#%>%
visIgraphLayout(layout = 'layout_nicely', smooth = TRUE)#%>%

#visLayout(randomSeed = 150 )
#improvedLayout =TRUE,
vis_net_vis
visSave(vis_net_vis, file = paste0(out_res,   '/vis_net.html'))





