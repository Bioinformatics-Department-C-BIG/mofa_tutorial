

## graph utils

#g<-OPI_g_de_mirs_de_genes_targets

get_logFC_by_node<-function(g){
    # get the colors for up and down for a graph 
    # @param
    # todo: add also an attribute size for the amount of DE? 
        g_names<-V(g)$name

        rna_fc<-rnas_sig[rnas_sig$GENE_SYMBOL %in% g_names , c('GENE_SYMBOL', 'log2FoldChange')]
        rna_fc$GENE_SYMBOL[rna_fc$log2FoldChange >0]
        mirs_fc<-mirnas_sig[mirnas_sig$GENE_SYMBOL %in% g_names , c('GENE_SYMBOL', 'log2FoldChange')]
        
        mirs_fc
        rna_fc
        fc_all=rbind(mirs_fc, rna_fc)

        up<-which(V(g)$name %in% mirs_fc$GENE_SYMBOL[mirs_fc$log2FoldChange >0])
        down<-which(V(g)$name %in% mirs_fc$GENE_SYMBOL[mirs_fc$log2FoldChange <0])

        g <- set.vertex.attribute(g, 'group',up , 'up')
        g <- set.vertex.attribute(g, 'group',down , 'down')

        up<-which(V(g)$name %in% rna_fc$GENE_SYMBOL[rna_fc$log2FoldChange >0]); up
        down<-which(V(g)$name %in% rna_fc$GENE_SYMBOL[rna_fc$log2FoldChange <0]); down

        g <- set.vertex.attribute(g, 'group',up , "up")
        g <- set.vertex.attribute(g, 'group',down , 'down')


        V(g)$color <- ifelse(igraph::V(g)$group == 'up', "#a02828","#126d8b")
        V(g)$FC<-fc_all[match(V(g)$name, fc_all$GENE_SYMBOL),]$log2FoldChange
        
        




    return(g)
}



library('visNetwork')


visualize_net<-function(visnet, net_name='net'){



  # visnet rectangle 
   visnet$nodes$font.size=35
   min(visnet$nodes$abs_FC*20, na.rm=TRUE)
   visnet$nodes$size= visnet$nodes$abs_FC*20

    visnet$nodes[is.na(visnet$nodes$abs_FC), ]$size = 5
    names(visnet$edges)
    visnet$edges<-visnet$edges[c('from', 'to')]
    vis_net_vis<-visNetwork(visnet$nodes, visnet$edges) %>%
               # visNodes( color =visnet$nodes$color  ) %>%
                visEdges(color='gray')

    vis_net_vis

    dir.create(paste0(outdir, '/networks/'))
    net_name=paste0('mirs_genes_', mofa_cluster_id, '_f',sel_factor,top_fr )
    net_name
    visSave(vis_net_vis, file = paste0(outdir, '/networks/',  net_name, '.html'))



    



    
   # vis_net_vis %>% visExport(type = "jpeg", name = "export-network", 
   #float = "left", label = "Save network", background = "white", style= "") 

   # visExport(  vis_net_vis,  file = "network.png", background = "black")

  return(vis_net_vis)
}


