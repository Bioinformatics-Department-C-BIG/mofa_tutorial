
source(paste0(script_dir,'ppmi/mofa_p_value_combination.R'))
source(paste0(script_dir,'ppmi/utils.R'))
source(paste0(script_dir,'multi_emap_mofa_plot.R'))

###########
# create multi omics emap plot 
library('multienrichjam')
library('clusterProfiler')
library('enrichplot')


which(f_pvals[[1]]$Description == 'regulation of small GTPase mediated signal transduction')
f_pvals[[1]][42,]



lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 


##############
########### MERGING



get_ranks<-function(f_pvals_merged_fs){
  #''
  #
  ranks_fs<-as.data.frame(apply(f_pvals_merged_fs, 2, rank))
  tot_rank=rowSums(ranks_fs)
  min_rank<-colnames(ranks_fs)[apply(ranks_fs,1,which.min)]
  min_rank
  rank_stats=data.frame(tot_rank=tot_rank, min_rank=min_rank)
  colnames(ranks_fs)<-paste0(colnames(ranks_fs), '_rank')
  
  rank_stats=cbind(rank_stats, ranks_fs)
  return(rank_stats)
}


f_pvals_merged<-f_pvals %>% purrr::reduce(full_join , by='Description' )
f_pvals_merged
dim(f_pvals_merged)

fish_cols<-colnames(f_pvals_merged)[ grepl('fish', colnames(f_pvals_merged )) & !grepl('pval', colnames(f_pvals_merged )) ]
f_pvals_merged_fs<-f_pvals_merged[,fish_cols]
colnames(f_pvals_merged_fs)<-names(sel_factors)
rownames(f_pvals_merged_fs)<-f_pvals_merged$Description
rownames(f_pvals_merged)<-f_pvals_merged$Description



### Apply Ranking: 


f_pvals_merged$Least_value<-rowMins(as.matrix(f_pvals_merged_fs))
f_pvals_merged$Least_factor<-colnames(f_pvals_merged_fs)[apply(f_pvals_merged_fs,1,which.min)]

# filter
f_pvals_merged_sig<-f_pvals_merged[f_pvals_merged$Least_value<0.05,]
f_pvals_merged_fs_sig<-f_pvals_merged_fs[f_pvals_merged$Least_value<0.05,]
### Get the ranks among the significant only !! 
rank_stats<-get_ranks(f_pvals_merged_fs_sig)
head(rank_stats)


#f_pvals_merged$tot_rank=rank_stats$tot_rank
#f_pvals_merged$min_rank=rank_stats$min_rank
f_pvals_merged_sig<-cbind(f_pvals_merged_sig,rank_stats)




########### FROM NOW ON USE SIG

rownames(f_pvals_merged_sig)<-standardize_go_names(rownames(f_pvals_merged_sig))
f_pvals_merged_sig$GOID<-go_ids$GOID[match(rownames(f_pvals_merged_sig),go_ids$TERM_standardized )]
f_pvals_merged_sig[, c('GOID','Least_factor')]


colnames(f_pvals_merged_sig)




rank_mofa_paths<-function(f_pvals){
        #' Rank the mofa pathways from the different factors. 
        #' TODO: weigh factors from each 
        #' @param merged_factors_mofa_sig description
        #' @param 
        #' TODO: redo
        #'
      
        ## what to actually input 
        # TODO: RANK factors? 
        lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 


        merge(f_pvals[[1]],f_pvals[[2]], by='Description' )
        f_pvals_merged<-f_pvals %>% reduce(full_join, by='Description')
        dim(f_pvals_merged)
        head(f_pvals_merged)
        
        
        f_pvals_merged$tot_rank=rowSums(f_pvals_merged[,c('rank1', 'rank2', 'rank3', 'rank4')])
        
        # weighted with vars
        ranks<-f_pvals_merged[,c('rank1', 'rank2', 'rank3', 'rank4')]
        vars_by_factor_sel<-rowSums(vars_by_factor)[sel_factors]
        vars_by_factor_sel<-vars_by_factor_sel/sum(vars_by_factor_sel)
        cors_by_factor_sel<-abs(cors_pearson_l[sel_factors,'CONCOHORT'])
                                
        f_pvals_merged$tot_rank
        
      
        ####        
        merged_factors_mofa<-do.call(rbind,f_pvals)
        #### Merge factors together??? #### 
        merged_factors_mofa_sig<-merged_factors_mofa[merged_factors_mofa$fish<0.05,]
        merged_factors_mofa_sig$Description<-gsub('-', ' ', tolower(merged_factors_mofa_sig$Description))
        unique(length(merged_factors_mofa_sig$Description))
  
        merged_results<-merged_factors_mofa_sig
        merged_results_o<-merged_results[
          with(merged_results, order(Description, fish,decreasing = FALSE)),
        ]
        merged_results_o[, c('fish', 'Description')]
        ## TODO:  double  IF i AM filteirng by largest
        merged_results<-merged_results_o[!duplicated(merged_results_o$Description),]
        merged_results[, c('fish', 'Description')]
        
        merged_results<-merged_results[
          with(merged_results, order(fish,decreasing = FALSE)),
        ]
        
        merged_results$pvalue<-merged_results$fish
        merged_results<-merged_results[merged_results$fish<0.05,]  

        
          
        
        
        
        return(list(merged_results, f_pvals_merged_ord))

        
        }

merged_results<-f_pvals_merged
write.csv(f_pvals_merged_sig, paste0(outdir, '/enrichment/merged_factors_pvals.csv'))

### TODO: where to obtain the gene names from??? 
# Now they are taken from RNA 
merged_results$geneColname<-enrich_rna[match(merged_results$Description, enrich_rna$Description),]$core_enrichment

enrich_res_df<-merged_results
colnames(merged_results)


############

create_enrich_result<-function(enrich_res_df){
  #' convert pvalue comb to enrich result for further cluster profiler plotting
  #' @param  enrich_res_df contains fish variable as pvalue
  #' @return description
  enr_full_all <- multienrichjam::enrichDF2enrichResult(as.data.frame(enrich_res_df),
                                                        keyColname =  'Description',
                                                        geneColname ='geneColname',
                                                        pvalueColname = 'fish',
                                                        descriptionColname = 'Description',
                                                        qvalue='',
                                                        pvalueCutoff = 0.05)
  return(enr_full_all)
}



N_EMAP=25
require(DOSE)
require(enrichplot)
require(viridis)
create_multi_emap<-function(enr_full, N_EMAP=25,fname_net, title, color_by='qvalue', min=0.2,
                            width=10, height=10, scale_option='plasma'){
          #'
          #'
          #'
          #'
          
          #N_EMAP=50
          #min=0.3
          #color_by='qvalue'
  
          options(ggrepel.max.overlaps = Inf)
          
          x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)
          
          
          
          
          
          
         # color_by='pvalue'
          #jpeg(fname_net)
          em<-emapplot(x2, showCategory=N_EMAP,
                       cex_category=1, 
                       cex_label_category=0.8, 
                       color=color_by,
                       min =min,
                       cluster.params=list(cluster =TRUE, 
                                method=cluster::kmeans, # cluster::clara, cluster::pam , cluster::kmeans
                                ellipse_style='ggforce'),
                      group_category=TRUE,
                    
                      alpha=0.1,
                      ellipse_style='ggforce'
                      # node_label = "all"
          )
          
          

          em
          
          em+ scale_fill_viridis(option=scale_option)
          #em +  scale_color_gradient(low = "#56B1F7", high = "#132B43")


          
          ggsave(fname_net, width=width, height=height, 
                 dpi = 600)  #ggsave(filename=fname_net, plot=last_plot(),width=15, height=15, dpi=300)
     
  
}






### Plot for each one independently 
use_mofa=TRUE
for (fn in fns){
        
        list_all<-get_mofa_paths_and_weights(factor=sel_factors[fn])
          # also write to extra file
          ### Create one emap for each factor with all paths inside factor 
        enrich_rna= list_all[[1]]
        enrich_proteins = list_all[[2]]
        enrich_mirnas = list_all[[3]]
        
        # outdir
          
        use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
        merged_path_file_mofa<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',run_ORA, pmethod,
                                        'mofa_',  use_mofa_s   )
        factor1_paths<-f_pvals[[fn]]
        factor1_paths=factor1_paths[factor1_paths$fish<0.05,]
        factor1_paths=factor1_paths[c('Description', 'fish')]
        #factor1_paths$rank=rank(factor1_paths$fish)
        ### remove low ranks ? 
        
        ## obtain the enrich RNA of the specific factor!!! 
        factor=names(sel_factors[fn])
        weights_Var=get_mofa_vars(factor=sel_factors[fn], adj_weights)
        
        
        factor1_paths$geneColname<-enrich_rna[match(factor1_paths$Description, enrich_rna$Description),]$core_enrichment
        max_ws<-names(which.max(weights_Var))
        max_ws

        
        enrich_res_df<-factor1_paths

        enr_full<-create_enrich_result(enrich_res_df)
        
        
        # run_enrichment_plots(enr_full,results_file = 'test.txt')
        enr_full@result$Description<-enr_full@result$ID
        
        
        min_overlap_score=0.25
        
        
        
    
        #### EMAP by factor ####
        create_multi_emap(enr_full, N_EMAP=50, fname_net=paste0(merged_path_file_mofa, '_network', max_ws, '.jpeg'), color_by = 'pvalue', 
                          min = min_overlap_score, scale_option='mako'
                            )
        
        print(merged_path_file_mofa)


}


#### Change the visualization ####

enr_full@result$Description<-enr_full@result$ID

library(arules)
library(igraph)
## create an emap as an igrpah object
N_EMAP=200
x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)

ig_plot<-enrichMapJam(x2, n=N_EMAP)
ig_plot
V(ig_plot)$size<-V(ig_plot)$size/1.2

V(ig_plot)$color<-factor(V(ig_plot)$min_rank)
V(ig_plot)$min_rank

#### Remove low rank ####
factor
high_quant_t=0.9
g_graph<-ig_plot





ig_plot_filt<-g_graph


ig_plot_filt<-remove_subcomponents(ig_plot_filt, subcomp_min_edge = 2)


#### Remove subcomponents 



plot.igraph(ig_plot_filt, 
     layout=layout_nicely, 
     
     
     # === vertex label
     vertex.label.family="Helvetica",                  # Font family of the label (e.g.“Times”, “Helvetica”)
     vertex.label.font=1,                         # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.7,                           # Font size (multiplication factor, device-dependent)
     
     )



### 1. Color by p-value, 
### 2.
vis_emap_mofa<-toVisNetworkData(ig_plot_filt)
vis_emap_mofa$edges$width
vis_emap_mofa$edges$length=1/vis_emap_mofa$edges$width*10



visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges) %>%
visIgraphLayout(layout = 'layout_nicely')# same as   visLayout(hierarchical = TRUE) 
  
head(vis_emap_mofa$edges$width)
vis_emap_mofa$nodes$id

visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges, 
           main=paste0(choose_f))%>%
  #visOptions(selectedBy= list(variable="group",multiple=T)) %>%
  #visEdges(width=width, 
  #         length=length)%>%
  visIgraphLayout(layout = 'layout_nicely') %>%# same as   visLayout(hierarchical = TRUE) 
  visNodes(font=list(size=20) )%>%
  visSave(file = paste0(outdir, '/enrichment/all_factors_emap_',pval_to_use, choose_f,'.html'))

#%>%



#### 2. Create one EMAP for all mofa factors ####
###


  list_all<-get_mofa_paths_and_weights(factor=sel_factors[fn])
  # also write to extra file
  # where to get lists of common genes 
  enrich_rna= list_all[[1]]
  enrich_proteins = list_all[[2]]
  enrich_mirnas = list_all[[3]]
  
  # outdir
  
  use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
  merged_path_file_mofa<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',run_ORA, pmethod,
                                'mofa_',  use_mofa_s   )
  
  merged_results=f_pvals_merged
  merged_results_sig=f_pvals_merged_sig
  ### ALL FACTOR PATHS 
  #merged_results_sig=merged_results[merged_results$Least_value<0.01,]
  dim(merged_results_sig)
  
  ### Create filters on nodes 
  merged_results_sig
  high_quant_t=0.5
  factor_filters<-names(sel_factors)
  rank_cols<-grep('[1-9]_rank', colnames(merged_results_sig))
  
  

  mark_high_rank<-function(x, high_quant_t){
    up_quant<-quantile(x, high_quant_t)
    return(x>up_quant)
  }
  
  
  merged_results_sig_fs<-merged_results_sig[,fish_cols]
  colnames(merged_results_sig_fs)<-names(sel_factors)
  
  rank_stats<-get_ranks(merged_results_sig_fs)
  merged_results_sig<-cbind(merged_results_sig,rank_stats)

  
  high_ranks<-apply(merged_results_sig[,rank_cols],2,mark_high_rank, high_quant_t=0.9)
  #high_ranks_3<-mark_high_rank(merged_results_sig$Factor3_rank,  high_quant_t=0.1)
  #high_ranks[,'Factor3_rank']<-high_ranks_3
  merged_results_sig<-merged_results_sig[!apply(high_ranks, 1,all),]
  dim(merged_results_sig)
  
  table(merged_results_sig$min_rank)
  
  #merged_results_sig=merged_results_sig[c('Description', 'fish')]
  ## obtain the enrich RNA of the specific factor!!! 

  ### Connect by gene OVERLAP 
  #merged_results_sig$Description
  #enrich_rna$Description
  merged_results_sig$geneColname<-enrich_rna[match(merged_results_sig$Description, enrich_rna$Description),]$core_enrichment
  merged_results_sig$Least_factor
  # hack to get it to color with this attribute 
  ## Make sure that pvalyus is the combined 
  merged_results_sig$pvalue= merged_results_sig$Least_value
  merged_results_sig$pvalue
  merged_results_sig$qvalue= merged_results_sig$min_rank
  merged_results_sig$qvalue=as.numeric(factor(merged_results_sig$qvalue))#, levels=c(1,2,3,4))
  enrich_res_df<-merged_results_sig

  enr_full<-create_enrich_result(enrich_res_df)
  enr_full@result$pvalue
  
  # run_enrichment_plots(enr_full,results_file = 'test.txt')
  enr_full@result$Description<-enr_full@result$ID
 # enr_full@result$qvalue<-as.factor(enr_full@result$qvalue)

  
    min_overlap_score=0.2
  #### EMAP by factor ####
  create_multi_emap(enr_full, N_EMAP=90, fname_net=paste0(merged_path_file_mofa, '_network_ALLfs', '.jpeg'),
                    color_by = 'qvalue', 
                    min = min_overlap_score, 
                    scale_option = 'plasma')
  
  
  















intersect(enr_full@result$Description, inter$a1)

### tidy up these stuff here: 
# show what is common in all
enr_full@result$qvalue[enr_full@result$Description %in% inter$a3]<-5 # rna only and mofa? 
enr_full@result$qvalue[enr_full@result$Description %in% inter$a7]<-1 # mofa unique
enr_full@result$qvalue[enr_full@result$Description %in% inter$a5]<-20 # all
enr_full@result$qvalue[enr_full@result$Description %in% inter$a5]<-20 # all
N_EMAP=100
graphics.off()
fname_net<-paste0(out_compare,'network_MOFA_single_',use_mofa_s, 
                  int_params ,N_EMAP,'.png')
create_multi_emap(enr_full, N_EMAP=50, fname_net=fname_net)
unique(inter$a7)



graphics.off()

create_multi_emap(enr_full, N_EMAP=100, fname_net=fname_net, min_edge = 0.99)

  




inter<-inter_mofa_union
unique(inter$a1)
#### Mofa vS single combination 
enr_full<-enr_full_all
enr_full@result$Description<-enr_full@result$ID
enr_full@result$qvalue[enr_full@result$Description %in% common_mofa_combination]<-10 # mirnas only

x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)

N_EMAP=50
jpeg(fname_net)
em<-emapplot(x2, showCategory=N_EMAP,
             cex_category=1, 
             cex_label_category=0.5, 
             color='qvalue'
)
em
dev.off()


create_multi_emap(enr_full, title='Mofa vs Union', fname_net = fname_net)
enr_full@result$qvalue

#### COMPARE SINGLE COMBINATION TO MOFA COMBINATION ####

cor_t<-0.15
outdir
#install.packages('arules')
library(arules)
library(igraph)
## create an emap as an igrpah object
N_EMAP=100
x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)

ig_plot<-enrichMapJam(x2, n=N_EMAP)

V(ig_plot)

grs=c(names(V(ig_plot)) %in% inter$a7)
grs[names(V(ig_plot)) %in% inter$a7]='mofa_only'
grs[names(V(ig_plot)) %in% inter$a5]='all'
grs[names(V(ig_plot)) %in% inter$a4]='rna_only'
grs[names(V(ig_plot)) %in% inter$a6]='mirna_only'

grs

names(V(ig_plot)) %in% inter$a7
ig_plot<-set_vertex_attr(ig_plot, 'group', index = V(ig_plot), grs)

library(RColorBrewer)

pal <- brewer.pal(length(unique(V(ig_plot)$group)), "Dark2")
pal
ig_plot
plot.igraph(ig_plot, 
     vertex.label.family='Helvetica', 
     vertex.color = pal[as.numeric(as.factor(vertex_attr(ig_plot, "group")))],
     )



my_object <- new("enrichResult",
                 readable = FALSE,
                 result = merged_results,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.2,
                 organism = "UNKNOWN",
                 ontology = "UNKNOWN",
                 gene = m,
                 keytype = "UNKNOWN",
                 universe = universe_vector,
                 gene2Symbol = character(0),
                 geneSets = geneSets)



