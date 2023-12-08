

## graph utils
g<-OPI_g_de_mirs_de_genes
V(g)$name

get_logFC_by_node<-function(g){
    # get the colors for up and down for a graph 
    # @param
    # todo: add also an attribute size for the amount of DE? 
        g_names<-V(g)$name

        rna_fc<-rnas_sig[rnas_sig$GENE_SYMBOL %in% g_names , c('GENE_SYMBOL', 'log2FoldChange')]
        rna_fc$GENE_SYMBOL[rna_fc$log2FoldChange >0]
        mirs_fc<-mirnas_sig[mirnas_sig$GENE_SYMBOL %in% g_names , c('GENE_SYMBOL', 'log2FoldChange')]
        up<-which(V(g)$name %in% mirs_fc$GENE_SYMBOL[mirs_fc$log2FoldChange >0])
        down<-which(V(g)$name %in% mirs_fc$GENE_SYMBOL[mirs_fc$log2FoldChange <0])

        g <- set.vertex.attribute(g, 'group',up , 'up')
        g <- set.vertex.attribute(g, 'group',down , 'down')

        up<-which(V(g)$name %in% rna_fc$GENE_SYMBOL[rna_fc$log2FoldChange >0]); up
        down<-which(V(g)$name %in% rna_fc$GENE_SYMBOL[rna_fc$log2FoldChange <0]); down

        g <- set.vertex.attribute(g, 'group',up , "up")
        g <- set.vertex.attribute(g, 'group',down , 'down')


        V(g)$color <- ifelse(igraph::V(g)$group == 'up', "pink","lightblue")
            


    return(g)
}

plot(g)

library('visNetwork')


visualize_net(toVisNetworkData(g))
V(g)
visnet=toVisNetworkData(g)
visnet$nodes

visualize_net<-function(visnet, net_name='net'){
    # visnet rectangle 
    visnet$nodes$font.size=30
    visnet$nodes$size=5
visnet$nodes$color
    vis_net_vis<-visNetwork(visnet$nodes, visnet$edges) %>%
                visNodes( color =visnet$nodes$color  ) %>%
                visEdges(color='gray')

    vis_net_vis

    dir.create(paste0(outdir, '/networks/'))
    visSave(vis_net_vis, file = paste0(outdir, '/networks/', net_name, '.html'))



    
   # vis_net_vis %>% visExport(type = "jpeg", name = "export-network", 
   #float = "left", label = "Save network", background = "white", style= "") 

   # visExport(  vis_net_vis,  file = "network.png", background = "black")

  return(vis_net_vis)
}
