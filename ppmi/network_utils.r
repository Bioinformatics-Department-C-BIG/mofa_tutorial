

## graph utils


rnas_V08$mofa_sign
#g<-OPI_g_de_mirs_de_genes_targets

get_logFC_by_node<-function(g, de_rnas, de_mirnas){
    # get the colors for up and down for a graph that is already built 
    # @param
    # todo: add also an attribute size for the amount of DE? 

    #
    #
    rnas_sig_visit<-de_rnas%>% dplyr::filter(mofa_sign == 'Significant')
    mirnas_sig_visit<-de_mirnas%>% dplyr::filter(mofa_sign == 'Significant')

        g_names<-V(g)$name

        # extract fold change from de results 
        rna_fc<-de_rnas[de_rnas$GENE_SYMBOL %in% g_names , c('GENE_SYMBOL', 'log2FoldChange')]
        mirs_fc<-de_mirnas[de_mirnas$GENE_SYMBOL %in% g_names , c('GENE_SYMBOL', 'log2FoldChange')]

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
        
        
        ## 

        V(g)$significant<-ifelse(V(g)$name %in% c(rnas_sig_visit$GENE_SYMBOL, 
                                  mirnas_sig_visit$GENE_SYMBOL), TRUE, FALSE)



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


  return(vis_net_vis)
}

library('OmnipathR')
create_regulatory_net_backbone<-function(rnas_sig_factor, mirnas_sig_factor, resources =c( "SIGNOR", "STRING_talklr" )){
  #'
  #' @param rnas_sig_factor significant de rnas
  #' @param mirnas_sig_factor
  #' take de rnas sig
  #' take de mirnas sign
  #' Load their interactions c
#
    ### Start loading dbs interactions, mirnas, mirtarges
        
      ## Until the DoRothEA issue gets fixed we have this here:
      interactions_string <-
          import_omnipath_interactions(resources=resources)


      # These will be used for the intermediates
      interactions_dor <- import_transcriptional_interactions(
          resources = resources
      )


      dim(interactions_dor)
      dim(interactions_string)

      ## ----mirnatarget--------------------------------------------------------------------------------------------
      ## We query and store the interactions into a dataframe
      interactions_mirs <-
        import_mirnatarget_interactions(resources = c("miR2Disease", "miRDeathDB", "miRTarBase", "TransmiR"))

      ## 1. interactions_mirs- obtain its genes 
      ## 2. filter by genes that have a de target OR are DE themselves 
      interactions_de_mirs<-interactions_mirs %>% 
                      dplyr::filter(source_genesymbol %in% mirnas_sig_factor$GENE_SYMBOL)

      dim(interactions_de_mirs)
      mirtargets<-interactions_de_mirs$target_genesymbol

      # find out which mir targets have a de interaction
      # TODO: this is only one way ie. if the mir is a source. 
      interactions_of_mirtargets_tar_is_de<-interactions_string %>% dplyr::filter(
          source_genesymbol %in% mirtargets) %>% dplyr::filter(
          target_genesymbol %in% rnas_sig_factor$GENE_SYMBOL

      )


      interactions_of_mirtargets_tar_is_de
      # 2. Now do the filter of gene targets 
      # 1. either de genes (in the network ) or de themselves 

      interactions_de_mirs_targets_filt<-interactions_de_mirs %>%
              dplyr::filter( target_genesymbol %in% interactions_of_mirtargets_tar_is_de$source_genesymbol | 
              target_genesymbol %in% rnas_sig_factor$GENE_SYMBOL # CAREFUL to access these with gene symbol!!
              )



      interactions_mirs$source_genesymbol
      interactions_de_mirs_targets_filt$source_genesymbol[grep( 'miR-7-5p',interactions_de_mirs_targets_filt$source_genesymbol)]
      interactions_de_mirs_targets_filt$target_genesymbol[grep( 'SNCA',interactions_de_mirs_targets_filt$target_genesymbol)]

      g_mirs_targets_filt<-interaction_graph(interactions_de_mirs_targets_filt)
      g_targets<-interaction_graph(interactions_of_mirtargets_tar_is_de)
      
      ############# Add also gene-gene interactions that are not in the graph already ###

      interactions_dor_target_genes<-interactions_dor %>%
          dplyr::filter(target_genesymbol %in% c(rnas_sig_factor$GENE_SYMBOL ) ) %>%
          dplyr::filter( source_genesymbol %in% c(rnas_sig_factor$GENE_SYMBOL )) 


      dim(interactions_dor_target_genes)
      g_genes_inter<-interaction_graph(interactions_dor_target_genes)


      g_extended<-union(g_mirs_targets_filt, g_targets) # mir-gene-gene interactions
      #g_extended<-union(g_extended, g_genes_inter) ## add gene-gene interactions

      return(g_extended)

}


library(multiMiR)

create_regulatory_net_backbone_mirtar<-function(rnas_sig_factor, mirnas_sig_factor){


    top_miRNAs<-mirnas_sig_factor$GENE_SYMBOL

  # Plug miRNA's into multiMiR and getting validated targets
  multimir_results <- get_multimir(org     = 'human',
                                 mirna   = top_miRNAs,
                                 table   = 'validated',
                                 summary = TRUE)

  head(multimir_results@data, 10)

  multimir_results_filt_targets<-multimir_results@summary[multimir_results@summary$target_symbol %in% 
                rnas_sig_factor_all_clusts$GENE_SYMBOL,]

  mirna_targets<-multimir_results_filt_targets[, c('mature_mirna_id','target_symbol' )]


  mirna_targets
  g<-graph_from_edgelist(as.matrix(mirna_targets), directed = TRUE)

return(g)
}

