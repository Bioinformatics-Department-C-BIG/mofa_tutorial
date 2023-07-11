
source(paste0(script_dir,'ppmi/mofa_p_value_combination.R')) ## run to merge 
source(paste0(script_dir,'ppmi/utils.R'))
source(paste0(script_dir,'multi_emap_mofa_plot.R'))
source(paste0(script_dir,'ppmi/emap_utils.R'))

###########
# create multi omics emap plot 
library('multienrichjam')
library('clusterProfiler')
library('enrichplot')
library('DOSE')



which(f_pvals[[1]]$Description == 'regulation of small GTPase mediated signal transduction')
f_pvals[[1]][42,]



lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 

##############
########### MERGING



get_ranks<-function(f_pvals_merged_fs){
  #''
  #' get the ranks of each factor, min_rank and total rank
  #' @param f_pvals_merged_fs column for each factor pvalue 
  #' @returns rank_stats
  #'
  ranks_fs<-as.data.frame(apply(f_pvals_merged_fs, 2, rank, ties.method='first'))
  tot_rank=rowSums(ranks_fs)
  min_rank<-colnames(ranks_fs)[apply(ranks_fs,1,which.min)]
  rank_stats=data.frame(tot_rank=tot_rank, min_rank=min_rank)
  colnames(ranks_fs)<-paste0(colnames(ranks_fs), '_rank')
  rank_stats=cbind(rank_stats, ranks_fs)
  return(rank_stats)
}


### Merge the pvalues for the factors  factors by pathway name 
f_pvals_merged<-f_pvals %>% purrr::reduce(full_join , by='Description' )
fish_cols<-colnames(f_pvals_merged)[ grepl('fish', colnames(f_pvals_merged )) & !grepl('pval', colnames(f_pvals_merged )) ]
f_pvals_merged_fs<-f_pvals_merged[,fish_cols]
colnames(f_pvals_merged_fs)<-names(sel_factors)
rownames(f_pvals_merged_fs)<-f_pvals_merged$Description
rownames(f_pvals_merged)<-f_pvals_merged$Description


pvalueT_mofa=0.01
### Apply Ranking: ####
f_pvals_merged$Least_value<-rowMins(as.matrix(f_pvals_merged_fs)) ## lower pvalue for path from all factors 
f_pvals_merged$Least_factor<-colnames(f_pvals_merged_fs)[apply(f_pvals_merged_fs,1,which.min)] 

# filter
f_pvals_merged_sig<-f_pvals_merged[f_pvals_merged$Least_value<pvalueT_mofa,]
f_pvals_merged_fs_sig<-f_pvals_merged_fs[f_pvals_merged$Least_value<pvalueT_mofa,]
### Get the ranks among the significant only !! 
rank_stats<-get_ranks(f_pvals_merged_fs_sig)
head(rank_stats)
### add the ranking statistics 
f_pvals_merged_sig<-cbind(f_pvals_merged_sig,rank_stats)




########### FROM NOW ON USE SIGnificant only ####

rownames(f_pvals_merged_sig)<-standardize_go_names(rownames(f_pvals_merged_sig))
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

merged_results<-f_pvals_merged ### THE MERGED mofa results 
write.csv(f_pvals_merged_sig, paste0(outdir, '/enrichment/merged_factors_pvals',pvalueT_mofa,  '.csv')) ## write the significant paths 

### TODO: where to obtain the gene names from??? 
# Now they are taken from RNA 
merged_results$geneColname<-enrich_rna[match(merged_results$Description, enrich_rna$Description),]$core_enrichment

enrich_res_df<-merged_results
colnames(merged_results)


############

create_enrich_result<-function(enrich_res_df, pvalueColname='fish', qvalue=''){
  #' convert pvalue comb to enrich result for further cluster profiler plotting
  #' @param  enrich_res_df contains fish variable as pvalue
  #' @return description
  enr_full_all <- multienrichjam::enrichDF2enrichResult(as.data.frame(enrich_res_df),
                                                        keyColname =  'Description',
                                                        geneColname ='geneColname',
                                                        pvalueColname = pvalueColname,
                                                        descriptionColname = 'Description',
                                                        qvalue=qvalue,
                                                        pvalueCutoff = 0.05)
  return(enr_full_all)
}



N_EMAP=25
require(DOSE)
require(enrichplot)
require(viridis)
library('stats')
library('cluster')
enr_full=enr_full_all_fs
create_multi_emap<-function(enr_full, N_EMAP=25,fname_net, title, color_by='qvalue', min=0.2,
                            width=10, height=10,
                            filter_single_nodes=FALSE, cex_label_category=0.6, nCluster=20,
                            node_label='category',scale_option='plasma', max.overlaps=10, 
                            layout='fr'){
          #'
          #' @param enr_full
          #' https://github.com/YuLab-SMU/enrichplot/blob/devel/R/emapplot.R
          #'
          #'
          
          #N_EMAP=2000;min=0.3;color_by='qvalue';cex_label_category=0.6; nCluster=50; filter_single_nodes=TRUE; max.overlaps=1
  
          options(ggrepel.max.overlaps =max.overlaps)
          #enr_full@result$ID[c(TRUE, TRUE,TRUE, FALSE)]<-'.'
          #enr_full@result$ID
          
          x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)
          x2_pass<-x2@termsim>min_overlap_score
          pos_edges<-rowSums(x2_pass, na.rm = TRUE) > 2L

          ### update without single nodes 
          enr_full_new<-enr_full[pos_edges, asis=T]
          #enr_full_new$Description=='organ growth'
          x2_new <- pairwise_termsim(enr_full_new, showCategory = N_EMAP)
          
          
         #grep( 'inflammatory',x2@result$ID )
         
         #enr_full@result[grep( 'inflammatory',enr_full@result$ID),]
         #enr_full@result$ID
          
          
          
          
          
         # color_by='pvalue'
          #jpeg(fname_net)
          
          if (filter_single_nodes){
            x2_to_use=x2_new
          }else{
            x2_to_use=x2
          }
        
          set.seed(150)
          em<-emapplot(x2_to_use, showCategory=N_EMAP,
                       cex_category=0.5, 
                       cex_label_category=cex_label_category,
                       color=color_by,
                       min =min,
                       cluster.params=list(cluster =TRUE, 
                                method=cluster::pam, # cluster::clara, cluster::pam , cluster::kmeans
                                ellipse_style='shadowtext', 
                                nCluster=nCluster, 
                               legend=TRUE),
                       layout.params =list(layout=layout),
                      group_category=TRUE,
                      nCluster=nCluster,
                      alpha=0.009,
                      ellipse_style='ggforce',
                      node_label = node_label, 
                      hilight.params=list(alpha_no_hilight =0.3, 
                                          alpha_hilight=0.3
                        
                      ),
          )
          
          em$data$name[!em$data$name %in% top_paths]<-''
          
          em$data$name[c(TRUE, TRUE, FALSE)]<-''
          
          
          em+ scale_fill_viridis(option=scale_option)
          #em +  scale_color_gradient(low = "#56B1F7", high = "#132B43")


          
          ggsave(fname_net, width=width, height=height, 
                 dpi = 600)  #ggsave(filename=fname_net, plot=last_plot(),width=15, height=15, dpi=300)
     
  
}






### Plot for each one independently 
use_mofa=TRUE
for (fn in fns){
        
        list_all<-get_mofa_paths_and_weights(factor=sel_factors[fn])
        print(sel_factors[fn])
          # also write to extra file
          ### Create one emap for each factor with all paths inside factor 
        enrich_rna= list_all[[1]]
        enrich_proteins = list_all[[2]]
        enrich_mirnas = list_all[[3]]
        
        # outdir
        #View(data.frame(enrich_mirnas$geneID))
        use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
        merged_path_file_mofa<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',run_ORA, pmethod,
                                        'mofa_',  use_mofa_s   )
        factor1_paths<-f_pvals[[fn]]
        
        ### leave only significant 
        factor1_paths=factor1_paths[factor1_paths$fish<0.05,]
        factor1_paths=factor1_paths[c('Description', 'fish')]
        #factor1_paths$rank=rank(factor1_paths$fish)
        ### remove low ranks ? 
        
        ## obtain the enrich RNA of the specific factor!!! 
        factor=names(sel_factors[fn])
        weights_Var=get_mofa_vars(factor=sel_factors[fn], adj_weights)
        
        
        factor1_paths$geneColname<-enrich_rna[match(factor1_paths$Description, enrich_rna$Description),]$core_enrichment
        max_ws<-names(which.max(weights_Var))
        min_overlap_score=0.25
        
        if (weights_Var['proteomics']>10){
          prot_names<-enrich_proteins[match(factor1_paths$Description, enrich_proteins$Description),]$core_enrichment
          length(prot_names)
          factor1_paths$geneColname<-enrich_rna[match(factor1_paths$Description, enrich_rna$Description),]$core_enrichment
          factor1_paths$proteinColname<-enrich_proteins[match(factor1_paths$Description, enrich_proteins$Description),]$core_enrichment
          factor1_paths$geneColname<-paste(factor1_paths$geneColname, factor1_paths$proteinColname, sep='/')
          min_overlap_score=0.25
          
          
        }else if(weights_Var['miRNA']>20){
          print(weights_Var)
          factor1_paths$geneColname<- enrich_mirnas[match(factor1_paths$Description, enrich_mirnas$Description),]$geneID
          min_overlap_score=0.60
          
        }
        
        
        ### TODo: Calculate a weighted overlap
        #else if (weights_Var['miRNA']>10){
        #  factor1_paths$
          #min_overlap_score=0.995
          
          
        #}


        enrich_res_df<-factor1_paths

      
        
        color_by='pvalue'
        color_by_rna=FALSE
        if (color_by_rna){
          
          factor1_paths$rna<-standardize_go_names(factor1_paths$Description) %in% standardize_go_names(mofa_rna_common)
        
          factor1_paths$qvalue=as.numeric(factor(factor1_paths$rna))
          factor1_paths[, c('rna', 'qvalue')]
        }
        enr_full<-create_enrich_result(enrich_res_df, pvalueColname = 'fish', qvalue =qvalue )
        
        # run_enrichment_plots(enr_full,results_file = 'test.txt')
        enr_full@result$Description<-enr_full@result$ID
        
        enr_full@result$Description
        factor1_paths$Description
        scale_option='mako'
        if (color_by_rna){

          enr_full@result$qvalue<-factor1_paths$qvalue
          color_by='qvalue'
          scale_option='plasma'
          
        }
        #### EMAP by factor ####
        create_multi_emap(enr_full, N_EMAP=80, fname_net=paste0(merged_path_file_mofa, '_network', max_ws,'rna' ,color_by_rna,'.jpeg'), 
                          color_by = color_by,  min = min_overlap_score, cex_label_category = 0.8,
                          nCluster=10, scale_option=scale_option, node_label = 'category', max.overlaps=20,
                          
                            )
       
        
        print(merged_path_file_mofa)


}








#### 2. Create one EMAP for all mofa factors ####
###

bl_f<-'ppmi/plots/p_BL_Plasma_0.9_T_1-2INEXPDvsn_TNA_0.9g_0.3_100_m_0.5_10_15_sig_FALSEcompleteFALSE_coh_1-2_BL_TRUE_split_FALSE/enrichment/merged_factors_pvals.csv'
bl<-read.csv(bl_f)
bl<-bl[!is.na(bl$Description),]
bl_sig<-bl[bl$Least_value<0.01,]

sel_factors[fn]

##TODO: merge the geneColnames from all factors..?? 
#https://rdrr.io/github/YuLab-SMU/enrichplot/src/R/emapplot_utilities.R

list_all<-get_mofa_paths_and_weights(factor=sel_factors[3])
list_all2<-get_mofa_paths_and_weights(factor=sel_factors[4])
sel_factors[2]
list_all3<-get_mofa_paths_and_weights(factor=sel_factors[1])
# also write to extra file
# where to get lists of common genes 
enrich_rna= list_all[[1]]
#enrich_proteins = list_all[[2]]
#enrich_mirnas = list_all[[3]]

# outdir

use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
merged_path_file_mofa<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',run_ORA, pmethod,
                              'mofa_',  use_mofa_s   )

merged_results=f_pvals_merged
merged_results_sig_all=f_pvals_merged_sig
### ALL FACTOR PATHS 
#merged_results_sig=merged_results[merged_results$Least_value<0.01,]
dim(merged_results_sig_all)

### Create filters on nodes 
high_quant_t=0.3
#factor_filters<-names(sel_factors)
rank_cols<-colnames(merged_results_sig_all)[grep('[1-9]_rank', colnames(merged_results_sig_all))]
rank_cols


mark_high_rank<-function(x, high_quant_t){
  up_quant<-quantile(x, high_quant_t)
  return(x>up_quant)
}


merged_results_sig_fs<-merged_results_sig_all[,fish_cols]
colnames(merged_results_sig_fs)<-names(sel_factors)

#rank_stats<-get_ranks(merged_results_sig_fs)
#cbind(merged_results_sig_fs$Factor1, rank(rank(merged_results_sig_fs$Factor1)))
# merged_results_sig<-cbind(merged_results_sig,rank_stats)
#dim(merged_results_sig)
table(merged_results_sig_all$min_rank)

### Keep top 20 paths from each?


## Filter more to show top representatives from each factor 
high_ranks<-apply(merged_results_sig_all[,rank_cols],2,mark_high_rank, high_quant_t=high_quant_t)
high_ranks
#high_ranks_3<-mark_high_rank(merged_results_sig$Factor3_rank,  high_quant_t=0.1)
#high_ranks[,'Factor3_rank']<-high_ranks_3
## remove low ranked in all 

merged_results_sig<-merged_results_sig_all[!apply(high_ranks, 1,all),]

dim(merged_results_sig)

## OR remove 50% in ALL factors ..?? 

### Add more genes? 
table(merged_results_sig$min_rank)
# 
factor_genes<-enrich_rna[match(merged_results_sig$Description, enrich_rna$Description),]$core_enrichment
factor_genes2<-list_all2[[1]][match(merged_results_sig$Description, list_all2[[1]]$Description),]$core_enrichment
factor_genes3<-list_all3[[1]][match(merged_results_sig$Description, list_all3[[1]]$Description),]$core_enrichment

ids_merge_all=c()
length(factor_genes)
for (i in length(factor_genes)){
  
  ids <- unique(unlist(strsplit(factor_genes[i], "/")))
  ids2 <-unique(unlist(strsplit(factor_genes2[i], "/")))
  ids3 <-unique(unlist(strsplit(factor_genes3[i], "/")))
  
  ids_merge<-union(ids,ids2)
  ids_merge
  ids_merge<-as.character(union(ids_merge,ids3))
  ids_merge
  ids_merge_all[i]<-paste0(ids_merge, sep = '/')
}

factor_genes_all<-paste0(factor_genes, factor_genes2, sep='/')
factor_genes_all<-paste0(factor_genes_all, factor_genes3, sep='/')
## filter to add only high var?? 
## Color genes by weight..?
merged_results_sig$geneColname<-factor_genes
merged_results_sig$geneColname<-factor_genes_all

# Color the nodes by qfactor
# hack to get it to color with this attribute 
## Make sure that pvalues is the combined --> here i use the least /minimum
# 1. color by factor where it has the minimum ranking 
#
merged_results_sig$pvalue = merged_results_sig$Least_value
merged_results_sig$qvalue = factor(merged_results_sig$min_rank)
merged_results_sig$qvalue = as.numeric(factor(merged_results_sig$qvalue))#, levels=c(1,2,3,4))

## OR color by whether it exists in baseline!!
color_by_rna=FALSE
color_by_bl=FALSE

min_overlap_score=0.25

if (color_by_bl){
  v08_in_BL<-intersect(bl_sig$Description, merged_results_sig$Description)
  merged_results_sig$BL<-merged_results_sig$Description %in% v08_in_BL
  merged_results_sig$qvalue=as.numeric(factor(merged_results_sig$BL))
  merged_results_sig[, c('BL', 'qvalue')]
}else if (color_by_rna){
  #  mofa_in_rna<-mofa_rna_common
  
  merged_results_sig$rna<-standardize_go_names(merged_results_sig$Description) %in% standardize_go_names(mofa_rna_common)
  merged_results_sig$qvalue=as.numeric(factor(merged_results_sig$rna))
  merged_results_sig[, c('rna', 'qvalue')]
  min_overlap_score=0.3
  
  
}

# TODO: color by miRNA unique or RNA unique

enrich_res_df<-merged_results_sig
enr_full<-create_enrich_result(enrich_res_df, pvalueColname = 'Least_value',qvalue = qvalue)
enr_full_all_fs<-enr_full
enr_full@result$pvalue


write.csv(merged_results_sig, paste0(merged_path_file_mofa, '_filtered_paths_network.csv'))

# run_enrichment_plots(enr_full,results_file = 'test.txt')
enr_full@result$Description<-enr_full@result$ID
# enr_full@result$qvalue<-as.factor(enr_full@result$qvalue)
top_paths<-enrich_res_df$Description[order(enrich_res_df$tot_rank, decreasing=FALSE) ][1:300]
length(enrich_res_df$Description)
top_paths

#### EMAP by factor ####
merged_path_file_mofa
node_label= 'group';cex_label_category=0.7
node_label= 'all';cex_label_category=0.5
node_label= 'category';cex_label_category=0.7

fname_net<-paste0(merged_path_file_mofa, '_network_ALLfs_', 'hq_', high_quant_t, 
                  'bl_',   as.numeric(color_by_bl)[1], 'rna_', as.numeric(color_by_rna)[1], 
                  node_label,'.jpeg')

 
create_multi_emap(enr_full, N_EMAP=1000, fname_net=fname_net,
                  width=15, height=10,
                  cex_label_category = cex_label_category,
                  filter_single_nodes=TRUE, 
                  nCluster=20,
                  max.overlaps=2,
                  ### TODO: add ggrob to add legend etc. 
                  color_by = 'qvalue', 
                  min = min_overlap_score,
                  #node_label = 'all',
                  node_label = node_label,
                  
                  scale_option = 'plasma')





####  TODO: FIND out which nodes to not have a single neighbour and delete them? 









#### Change the visualization ####

enr_full@result$Description<-enr_full@result$ID

library(arules)
library(igraph)
## create an emap as an igrpah object
N_EMAP=100
x2 <- pairwise_termsim(enr_full_all_fs, showCategory = N_EMAP)
enr_full_all_fs
ig_plot<-enrichMapJam(x2, n=N_EMAP)
# set overlap score
V(ig_plot)$size<-V(ig_plot)$size/1.8

V(ig_plot)$min_rank

vis_emap_mofa<-toVisNetworkData(ig_plot)
vis_emap_mofa$edges$width

vis_emap_mofa$edges$weight=width
vis_emap_mofa$edges$length=1/vis_emap_mofa$edges$width*10



visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges) %>%
  visLayout(hierarchical = TRUE) 

#visIgraphLayout(layout = 'hierarchical')# same as   visLayout(hierarchical = TRUE) 





#### Remove low rank ####
factor
high_quant_t=0.9
g_graph<-ig_plot
ig_plot_filt<-g_graph
#### Remove subcomponents 
hist(E(ig_plot_filt)$overlap)
which(E(ig_plot_filt)$overlap <0.25)


clust_labs<-NULL
### MIN overlap >0.35
ig_plot_filt <- delete.edges(ig_plot_filt, which(E(ig_plot_filt)$overlap <0.3) )
ig_plot_filt


filt_plot<-remove_subcomponents(ig_plot_filt, subcomp_min_edge = 1)
ig_plot_filt=filt_plot[[1]]

cluster_louvain(ig_plot_filt)
### Clustering

#clust_labs<-igraph::cluster_optimal(ig_plot_filt )
clust_labs_low_res<-cluster_louvain(ig_plot_filt, resolution=0.01)
clust_labs<-cluster_louvain(ig_plot_filt, resolution=1)

LO = layout_with_fr(ig_plot_filt)
## identify which communities have fewer than 5 members

members<-melt(clust_labs$memberships)
colnames(members)=c('ni','id','group')

group_ids <- lapply(members %>% split(.$group), function(grp) { grp$id })
group_ids
group_ids
ig_plot_filt
main_labels<-lapply(group_ids, function(x){if (length(x)>3)
                                            return(c(max(x), min(x))) })
main_label=unlist(main_labels, use.names=FALSE)
#main_label=unlist(main_label, use.names=FALSE)
label=V(ig_plot_filt)$label


#### Remove some of the labels ########### 
label=V(ig_plot_filt)$label
label
ig_plot_filt_2=ig_plot_filt
V(ig_plot_filt)$label<-V(ig_plot_filt)$name
label[!seq(1:NROW(label)) %in% main_label]<-''
label


cols_pal<-c("#bef7ff", "#a0dcff", "#82c2ff", "#63a7ff")
cols_pal<-RColorBrewer::brewer.pal(4, name='Spectral') # red, 2: orange 3: green , blue

fnames<-unique(factor(V(ig_plot_filt_2)$min_rank))

mapdf <- data.frame(old=fnames,new=cols_pal[1:length(sel_factors)])
V(ig_plot_filt_2)$color<-mapdf$new[match(V(ig_plot_filt_2)$min_rank, mapdf$old)]
V(ig_plot_filt_2)$color
ig_plot_filt_2 <- set.vertex.attribute(ig_plot_filt_2, "label", value=label)
#ig_plot_filt_2 <- set.vertex.attribute(ig_plot_filt_2, "color", value=)

V(ig_plot_filt_2)$label=label

choose_f=14
svg(paste0(outdir, '/enrichment/all_factors_emap_',pval_to_use, choose_f,'.svg'))
plot.igraph(ig_plot_filt_2, 
     layout=layout_nicely ,
     mark.groups = clust_labs,
     mark.col=c("tan", "tan","pink", "lightgray"), 
     
     # === vertex label
     vertex.label.family="Helvetica",                  # Font family of the label (e.g.“Times”, “Helvetica”)
     vertex.label.font=1,                         # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.cex=0.8                          # Font size (multiplication factor, device-dependent)    )
)
dev.off()



### 1. Color by p-value, 
### 2.
vis_emap_mofa<-toVisNetworkData(ig_plot)
vis_emap_mofa$edges$width

vis_emap_mofa$edges$weight=width
vis_emap_mofa$edges$length=1/vis_emap_mofa$edges$width*10

visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges) %>%
  visLayout(hierarchical = TRUE) 

      visIgraphLayout(layout = 'hierarchical')# same as   visLayout(hierarchical = TRUE) 
  
head(vis_emap_mofa$edges$width)
vis_emap_mofa$nodes$id

choose_f=sel_factors[fn]
visNetwork(vis_emap_mofa$nodes, vis_emap_mofa$edges, 
           main=paste0(choose_f))%>%
  #visOptions(selectedBy= list(variable="group",multiple=T)) %>%
  #visEdges(width=width, 
  #         length=length)%>%
  visIgraphLayout(layout = 'layout_nicely') %>%# same as   visLayout(hierarchical = TRUE) 
  visNodes(font=list(size=20) )%>%
  visSave(file = paste0(outdir, '/enrichment/all_factors_emap_',pval_to_use, choose_f,'.html'))

#%>%












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
x2 <- pairwise_termsim(enr_full_all_fs, showCategory = N_EMAP)

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



