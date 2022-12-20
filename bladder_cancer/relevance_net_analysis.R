
library(igraph)
cor_net<-correlation_net$gR
cor_net
is_bipartite(cor_net)

## TODO: negative and positive weights--> how to plot? just make them different color/attribute
plot(cor_net, 
     arrow.mode=0,
     vertex.size=20)


## Make a projection to monopartite