### Load the rnas and mirnas

# 1. 
#BiocManager::install('OmnipathR')
#BiocManager::install('multiMiR')

library('OmnipathR')

library(OmnipathR)
library(tidyr)
#library(dnet)
#BiocManager::install('dnet')
#install.packages('dnet')
#library(gprofiler2)
# interactions - proteins
library('igraph')
mofa_cluster_id<-2
# TODO: fix to choose the visit to plot.... 
# and the backbone visit - should be latest and all other to use that 
VISIT_COMP<-'V08'
padj_T<-0.5
cell_corr_deseq=FALSE
print(outdir )
#k_centers_m=3; rescale_option=TRUE; sel_group_cors = FALSE
#remove_cell_factors=FALSE

fact<-get_factors_for_metric(DIFF_VAR)
fact_s<-paste0(fact, collapse='_')
cluster_params_to_load<-cluster_params_dir
cluster_params_dir


load_de_by_visit_and_cluster<-function(VISIT_COMP , mofa_cluster_id, cell_corr_deseq = TRUE, prefix='rnas_' ){
     #'
     #' @param VISIT_COMP
     #' @param mofa_cluster_id   
    #' @param cell_corr   

    rnas_visit<-read.csv(paste0(cluster_params_dir, '/de_c', as.numeric(cell_corr_deseq), '/' , VISIT_COMP ,  '/', prefix,  'de_cluster_', mofa_cluster_id, '.csv'))

    if (prefix == 'mirnas_'){
      rnas_visit$GENE_SYMBOL<-rnas_visit$X 
    }

    return(rnas_visit )

}
output_files



load_proteins_by_visit<-function(cluster_id=1,VISIT_COMP, prefix='prot_'){
  
  proteins_visit_f<-paste0(cluster_params_dir, '/de_c0/',VISIT_COMP, '/' ,
                                         prefix, tissue,'_', prot_de_mode,'_de_cl',cluster_id,  '_results.csv')
                proteins_visit<-read.csv(proteins_visit_f)


}



# proteins_V08




clust_ids<-c(1,2,3, '1_2_3')
clust_ids<-c('1')
clust_ids<-c(1,2,3,'1_2_3')

# $adj.P.Val





rnas_V08_list<-lapply(clust_ids, load_de_by_visit_and_cluster, VISIT_COMP='V08', cell_corr_deseq=cell_corr_deseq, prefix= 'rnas_') # holds all clusters 
rnas_sig_V08_list<-lapply(rnas_V08_list, function(df){ return(df%>% dplyr::filter(padj<padj_T))})
lapply(rnas_sig_V08_list, dim)

names(rnas_V08_list)<-clust_ids
names(rnas_sig_V08_list)<-clust_ids
#proteins_V08_list<-lapply(clust_ids,load_proteins_by_visit, VISIT_COMP='V08')

rnas_sig_V08_list

load_all_times=TRUE
if (load_all_times){
    # TODO: for loop
  rnas_V06_list<-lapply(clust_ids, load_de_by_visit_and_cluster, VISIT_COMP='V06', cell_corr_deseq=cell_corr_deseq, prefix= 'rnas_') # holds all clusters 
  rnas_sig_V06_list<-lapply(rnas_V06_list, function(df){ return(df%>% dplyr::filter(padj<padj_T))})
  names(rnas_V06_list)<-clust_ids

  mirnas_V06_list<-lapply(clust_ids, load_de_by_visit_and_cluster, VISIT_COMP='V06', cell_corr_deseq=cell_corr_deseq, prefix= 'mirnas_')
  mirnas_sig_V06_list<-lapply(mirnas_V06_list, function(df){  return(df%>% dplyr::filter(padj<padj_T))})
  names(mirnas_V06_list)<-clust_ids


  rnas_BL_list<-lapply(clust_ids, load_de_by_visit_and_cluster, VISIT_COMP='BL', cell_corr_deseq=cell_corr_deseq, prefix= 'rnas_') # holds all clusters 
  rnas_sig_BL_list<-lapply(rnas_BL_list, function(df){ return(df%>% dplyr::filter(padj<padj_T))})
  names(rnas_BL_list)<-clust_ids

  mirnas_BL_list<-lapply(clust_ids, load_de_by_visit_and_cluster, VISIT_COMP='BL', cell_corr_deseq=cell_corr_deseq, prefix= 'mirnas_')
  mirnas_sig_BL_list<-lapply(mirnas_BL_list, function(df){  return(df%>% dplyr::filter(padj<padj_T))})
  names(mirnas_BL_list)<-clust_ids



}


# TODO:  union of all visits 




## load mirs
mirnas_V08_list<-lapply(clust_ids, load_de_by_visit_and_cluster, VISIT_COMP=VISIT_COMP, cell_corr_deseq=cell_corr_deseq, prefix= 'mirnas_')
mirnas_sig_V08_list<-lapply(mirnas_V08_list, function(df){  return(df%>% dplyr::filter(padj<padj_T))})

names(mirnas_V08_list)<-clust_ids
names(mirnas_sig_V08_list)<-clust_ids
mirnas_sig_V08_list






#### CHOOSE VISIT and cluster HERE
mirnas_visit<-mirnas_V08_list[[mofa_cluster_id]] # chosen visit to analyse
mirnas_sig_visit<-mirnas_sig_V08_list[[mofa_cluster_id]] # chosen visit to analyse

### define the rnas fcs that will be used
 
rnas_visit <- rnas_V08_list[[mofa_cluster_id]]
rnas_sig_visit<-rnas_sig_V08_list[[mofa_cluster_id]]


##### filter mofa top


rnas_top<-rnas_visit %>%arrange(padj) %>%  dplyr::filter(padj<padj_T) %>% as.data.frame
mirnas_top<-mirnas_visit %>%arrange(padj) %>% dplyr::filter(padj<padj_T) %>% as.data.frame
mirnas_top
fact
#sel_factor=fact[4]; 
sel_factor=fact[4]; 

if (cell_corr_deseq){
  top_fr=0.15 # FACTOR 11
  top_fr=0.11
  
  }else{
        top_fr=0.05
  }
top_genes_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'RNA',  factors=sel_factor, top_fr = top_fr))
top_genes_factor8<-get_symbols_vector(top_genes_factor8)
top_fr_mirs=0.1
top_mirnas_factor8<-gsub('\\..*','',select_top_bottom_perc(MOFAobject, 'miRNA',  factors=sel_factor, top_fr = top_fr_mirs))

length(top_genes_factor8)

## 1. Colour if in TOP 1% of factor
## 2. Colour if in specific pathways 
## 
dim(rnas_top); dim(mirnas_top)
top_g<-10000
top_m<-200

top_rnas<-rnas_top[1:top_g,]
top_mirnas<-mirnas_top[1:top_m,]
top_mirnas$GENE_SYMBOL<-top_mirnas$X






source(paste0(script_dir, 'ppmi/network_utils.R'))
############# MY ADDITION 
# 1. DE MIRS 
# 2. DE GENES 
mirnas_sig_visit

mirnas_sig_factor<-mirnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8)
rnas_sig_factor<-rnas_sig_visit %>% dplyr::filter(GENE_SYMBOL %in% top_genes_factor8)

rnas_sig_factor
## compare clust 1,2


mirnas_sig_factor_V08_list<-lapply(mirnas_sig_V08_list, function(df){
  return( df %>%  dplyr::filter(GENE_SYMBOL %in% top_mirnas_factor8))
})

rnas_sig_factor_V08_list<-lapply(rnas_sig_V08_list, function(df){
  return( df %>%  dplyr::filter(GENE_SYMBOL %in% top_genes_factor8))
})

## input: 
## de_rnas: union of all
## de_mirnas: union of all clusters 
mofa_filter=TRUE; 
if (mofa_filter){
  rnas_sig_list = rnas_sig_factor_V08_list
  mirnas_sig_list = mirnas_sig_factor_V08_list
}else{
 rnas_sig_list = rnas_sig_V08_list
  mirnas_sig_list = mirnas_sig_V08_list
  top_fr=1*2 # helps in plotting...

}

names(rnas_sig_list)

rnas_sig_factor_all_clusts<-Reduce(function(x, y) merge(x, y, all=TRUE), rnas_sig_list)
mirnas_sig_factor_all_clusts<-Reduce(function(x, y) merge(x, y, all=TRUE), mirnas_sig_list)

mirnas_sig_factor_all_clusts$GENE_SYMBOL
resources = c(  "STRING_talklr", 'SIGNOR')
 gene_mir_resources = c("miRTarBase", "TransmiR")
add_gg_interactions = TRUE


g_extended_both_clusters<-create_regulatory_net_backbone(rnas_sig_factor_all_clusts, mirnas_sig_factor_all_clusts, resources = resources, add_gg_interactions=add_gg_interactions)

length(g_extended_both_clusters)

#colnames(proteins_V08_list[[1]])[2] = 'log2FoldChange' # rename to align with rna columns
#proteins_V08_list[[1]]$GENE_SYMBOL=proteins_V08_list[[1]]$X
#proteins_V08_list[[1]]$mofa_sign = proteins_V08_list[[1]]$sign_lfc
# choose the backbone to be the union or intersection
rnas_sig_list[[1]]$type='gene'
#proteins_V08_list[[1]]$type='protein'
colnames(proteins_V08_list[[1]])
#proteins_V08_list[[1]]
#proteins_V08_list[[1]]$mofa_sign<-proteins_V08_list[[1]]$adj.P.Val <0.05
#proteins_V08_list[[1]$mofa_sign



rnas_proteins_sig_list<-rbind(rnas_sig_list[[1]][c('log2FoldChange', 'GENE_SYMBOL', 'mofa_sign', 'type')], proteins_V08_list[[1]][c('log2FoldChange', 'GENE_SYMBOL', 'mofa_sign', 'type')] )


 # TODO: move OMNIPATH query to file earlier do not rerun 
g_extended_both_clusters<-create_regulatory_net_backbone(rnas_sig_list[[1]], mirnas_sig_list[[1]],resources = resources, 
                                                        gene_mir_resources =  gene_mir_resources ,
                                                        add_gg_interactions=add_gg_interactions)



g_extended_both_clusters<-create_regulatory_net_backbone(rnas_sig_factor_all_clusts, mirnas_sig_factor_all_clusts,resources = resources, 
                                                        gene_mir_resources =  gene_mir_resources ,
                                                        add_gg_interactions=add_gg_interactions)


rnas_sig_factor_all_clusts

### Now plot! 
g_extended_both_clusters
OPI_g_union<-g_extended_both_clusters
OPI_g_union
# TODO: give the layout 

set.seed(123) 
# color by logFC 

# backbone should be the union 
#g_fc_V08<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08,  de_mirnas=mirnas_V08_list[[3]])




# TODO: apply take mis and genes 


# Unified plot: 
g_fc_V08_list<-sapply(clust_ids, function(mofa_cluster_id){
  g_fc_V08_cl<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V08_list[[mofa_cluster_id]],  de_mirnas=mirnas_V08_list[[mofa_cluster_id]])
  return(g_fc_V08_cl)


})

names(g_fc_V08_list)<-clust_ids

g_fc_V08_list


if (load_all_times){
  g_fc_V06_list<-sapply(clust_ids, function(mofa_cluster_id){
  g_fc_V06_cl<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_V06_list[[mofa_cluster_id]],  de_mirnas=mirnas_V06_list[[mofa_cluster_id]])
  return(g_fc_V06_cl)
    })
# load BL 
  g_fc_BL_list<-sapply(clust_ids, function(mofa_cluster_id){
  g_fc_BL_cl<-get_logFC_by_node(OPI_g_union, de_rnas=rnas_BL_list[[mofa_cluster_id]],  de_mirnas=mirnas_BL_list[[mofa_cluster_id]])
  return(g_fc_BL_cl)
})



}


OPI_g_union




## igraph plotting 

### Plot with specified coords 
# Get the coordinates of the Nodes of the backbone that includes all clusters 

Coords <- layout_with_fr(OPI_g_union)# %>% 
  Coords
  #as_tibble %>%
  #  bind_cols(data_frame(names = names(V(g_fc))))



## Plot cluster 3 as specified by backbone coords 

dir.create(paste0(outdir,cluster_params_to_load, '/nets/mf_', as.numeric(mofa_filter), '/'), recursive = TRUE)

mofa_cluster_id=1
set.seed(123) 

# Plot each cluster on the backbone 


      #V(g_fc_plot)$size<-V(g_fc_plot)$size*
      plot_settings_omnipath<-list(
        label.cex=1,
        vertex.size.factor=1
      )
      plot_settings_mirtar= list(
        label.cex=1,
        vertex.size.factor=1 
      )



# TODO: choose backbone to be from the intersection or run with all 3...? 


height_offset = 4
width_offset = 3
hw_multiplier = 0.8
vertex.frame.width = 1 # size of outline that marks significance 
vertex.label.dist = 1.5

dir.create(paste0(cluster_params_to_load, '/de_c',as.numeric(cell_corr_deseq), '/nets/', 'mf_', as.numeric(mofa_filter) ,'/' ), recursive=TRUE)

for (mofa_cluster_id in clust_ids){


 
     g_fc_plot<-g_fc_V06_list[[mofa_cluster_id]];  VISIT_PLOT = 'V06'
     g_fc_plot<-g_fc_BL_list[[mofa_cluster_id]];  VISIT_PLOT = 'BL'
      g_fc_plot<-g_fc_V08_list[[mofa_cluster_id]];VISIT_PLOT = 'V08'

      ## ADD border if significant 
      # create a vertor of border colours conditional on node type

      #V(g_fc_V06_list[[mofa_cluster_id]])$FC

 
      #round(cbind(V(g_fc_V08_list[[mofa_cluster_id]])$FC,V(g_fc_V06_list[[mofa_cluster_id]])$FC), digits=2)

      sig_names<-V(g_fc_plot)$name[V(g_fc_plot)$significant]
      bd <- ifelse(V(g_fc_plot)$significant, "#2fc729", NA) 
      V(g_fc_plot)$size<-ifelse(!is.na(V(g_fc_plot)$FC), log2(abs(V(g_fc_plot)$FC)+1)*17, 2)
      g_fc_plot$layout<-Coords

      plot_settings=plot_settings_omnipath

      V(g_fc_plot)$size<-V(g_fc_plot)$size*plot_settings$vertex.size.factor

      ### IGRAPH TO SAVE 
      net_name=paste0('mirs_genes_',VISIT_PLOT, '_f',sel_factor,top_fr, 'gg_', as.numeric(add_gg_interactions ),'cl_',mofa_cluster_id)
      net_file=paste0(paste0(cluster_params_to_load, '/de_c',as.numeric(cell_corr_deseq), '/nets/mf_', as.numeric(mofa_filter)) ,'/' , net_name)
      V(g_fc_plot)$frame.color<-bd
      V(g_fc_plot)$label.cex=plot_settings$label.cex # increase label size? 


      graphics.off()


      height=height_offset+log2(length(g_fc_plot))*hw_multiplier;
      width=width_offset+log2(length(g_fc_plot))*hw_multiplier; # Set size based on 
      png(paste0(net_file, '.png'), res=200, units='in', height=height, width=width)
      plot(g_fc_plot, 
      vertex.color=(adjustcolor(V(g_fc_plot)$color, alpha.f = 0.8 )), 
      vertex.size=V(g_fc_plot)$size,
     # vertex.shape=ifelse(grepl('miR|hsa-let-',V(g_fc_plot)$name ), 'rectangle', 'circle' ),
     # vertex.shape=ifelse(V(g_fc_plot)$type=='protein' , 'rectangle', 'circle' ),

      vertex.frame.color=bd,
      vertex.frame.width=vertex.frame.width,
      vertex.label.dist =vertex.label.dist, 
        layout = Coords)

      title(paste(  'Resources: ', paste(as.character(resources), collapse=', '),'\n', 
       paste(as.character(gene_mir_resources), collapse=', '),'\n', 
      'gene-gene ints:', add_gg_interactions, 
      ', Visit: ', VISIT_PLOT, ', Cluster: ', mofa_cluster_id, 
      ', \nMOFA filter', mofa_filter, ', Factor: ',
      sel_factor,', top: ', top_fr  ),cex.main=1)

      #ggsave(paste0(net_file, '.png'))
      dev.off()



      #### VISNWT TO SAVE


      #V(g_fc)$shape<- ifelse(grepl("miR|hsa-let",igraph::V(g)$name), "vrectangle", "circle")
      visnet <- toVisNetworkData(g_fc_plot)
      visnet$nodes$abs_FC<- abs(visnet$nodes$FC)
      
     # visSave(visIgraph(g_fc_plot), file =paste0(net_file, '.html'))




}






all_sig<-do.call(cbind,lapply(g_fc_V08_list,function(g_fc_plot){as.numeric(V(g_fc_plot)$significant) }))

all_sig<-as.data.frame(all_sig)
rownames(all_sig)<-V(g_fc_plot)$name
all_sig
write.csv(all_sig,paste0(net_file, '.csv'))