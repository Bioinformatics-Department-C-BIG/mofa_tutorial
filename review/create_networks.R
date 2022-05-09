



source('utils.R')
source('literature_review.R')
install.packages('visNetwork')

library('visNetwork')
##### 
#' Create an edge list from the frequency table ! 
# aggregated or separately? 

## todo label the size of the nodes df_to_plotby frequency 
cancer_filter='yes'
comb_freq<- comb_frequencies_by_group %>% filter(Cancer ==cancer_filter)
single_omics_frequencies_filtered<-single_omics_frequencies %>% filter(Cancer ==cancer_filter)


# Add weight as the value of the combination frequency
edge_list<-data.frame(do.call(rbind, str_split(comb_freq$Var1, ' - ')))
edge_list$weight<-comb_freq$Freq
#edge_list<-edge_list[edge_list$weight>10, ] 
edge_list<-edge_list[order(edge_list$weight, decreasing = TRUE),]

# most frequent datasets
# most_frequent<-comb_freq[order(comb_freq$perc, 
#                                                decreasing = TRUE),]
# most_frequent_5<-most_frequent$Var1[1:5]
# edge_list

# Also add the value on the node 
write.csv(edge_list,paste0('edge_list_',cancer_filter,'.csv'), quote=FALSE)

library(igraph)

#### Add the frequency of single omics to the network vertices
df<-single_omics_frequencies_filtered

net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)
net_att<-df[match(V(net)$name, df$Var1),]
df$names<-
  cancer_filter
#V(net)$color='lightblue'
#V(net)$frame.color='lightblue'
V(net)$size=log2(net_att$Freq)*6
E(net)$width=log2(edge_list$weight)
E(net)$label= edge_list$weight
data <- toVisNetworkData(net)
data$nodes$color='lightblue'
data$nodes$size=log2(net_att$Freq)*3

data$edges$color='blue'

svg(filename = 'network.svg', width=20, height=20)

visNetwork(nodes = data$nodes, edges = data$edges) %>%
  visNodes(
    color = list(border='blue'),
    font='10',
    x=c(1,2,3,4,5,6), y=c(1,1,1,2,2,2))  %>%
    visEdges(label = as.character(edge_list$weight) ) %>%
visIgraphLayout(layout='layout_as_tree')

dev.off()

  
visIgraph(net) %>% 
  visIgraphLayout(layout='layout_as_tree')





# TODO: assign omics frequencies of all samples or only cancer? 
#vertex_attr(net, 'freq', index=V(net))<-single_omics_frequencies$Freq
g <- make_lattice( c(3,3) )
layout_on_grid(g)

svg(width=20, height=20)

p<-plot.igraph(net, edge.width=log2(edge_list$weight), 
               vertex.size=log2(net_att$Freq)*6, 
               edge.label= edge_list$weight, color.edge='blue', 
               layout=layout_as_tree, 
               vertex.frame.color='lightblue', 
               vertex.frame.width=10)



#vertex.shape='crectangle')
dev.off()

visIgraph(net)
#save(p,paste0('plots/network', as.character(colname), '.png'))


#g <- set.vertex.attribute(g,'id',1,'first_id')


