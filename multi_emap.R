

###########
# create multi omics emap plot 
library('multienrichjam')
library('clusterProfiler')
library('enrichplot')


colnames(merged_factors_mofa)
colnames(enrich_rna)
merged_factors_mofa[[1]]
sel_factors
merged_factors_mofa_sig
which(f_pvals[[1]]$Description == 'regulation of small GTPase mediated signal transduction')
f_pvals[[1]][42,]





rank_mofa_paths<-function(f_pvals){
        #' Rank the mofa pathways from the different factors. 
        #' TODO: weigh factors from each 
        #' @param merged_factors_mofa_sig description
        #' @param 
        #'
        #'
      
        ## what to actually input 
        # TODO: RANK factors? 
        lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 

        lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 
        
        merge(f_pvals[[1]],f_pvals[[2]], by='Description' )
        f_pvals_merged<-f_pvals %>% reduce(inner_join, by='Description')
        
        
        f_pvals_merged$rank1[order(f_pvals_merged$fish.x)]<-1:nrow(f_pvals_merged)
        f_pvals_merged$rank2[order(f_pvals_merged$fish.y)]<-1:nrow(f_pvals_merged)
        f_pvals_merged$rank3[order(f_pvals_merged$fish.x.x)]<-1:nrow(f_pvals_merged)
        f_pvals_merged$rank4[order(f_pvals_merged$fish.y.y)]<-1:nrow(f_pvals_merged)
        
        f_pvals_merged$tot_rank=rowSums(f_pvals_merged[,c('rank1', 'rank2', 'rank3', 'rank4')])
        
        f_pvals_merged[order(f_pvals_merged$tot_rank),'Description']
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

        
          
        
        
        
        return(merged_results)
}




### TODO: where to obtain the gene names from??? 
# Now they are taken from RNA 
merged_results$geneColname<-enrich_rna[match(merged_results$Description, enrich_rna$Description),]$core_enrichment
enr_full <- multienrichjam::enrichDF2enrichResult(as.data.frame(merged_results),
                                                  keyColname =  'Description',
                                                  geneColname ='geneColname',
                                                  pvalueColname = 'fish',
                                                  descriptionColname = 'Description',
                                                  qvalue='',
                                                  pvalueCutoff = 0.05)

enr_full@result$pvalue



inter<-inter_mofa_single
# run_enrichment_plots(enr_full,results_file = 'test.txt')
enr_full@result$ID
enr_full@result$Description<-enr_full@result$ID
#enr_full@result$qvalue[enr_full@result$Description %in% inter$a6]<-10 # mirnas only
enr_full@result$qvalue[enr_full@result$Description %in% inter$a4]<-5 # rna only and mofa? 
enr_full@result$qvalue[enr_full@result$Description %in% inter$a7]<-1 # mofa unique
enr_full@result$qvalue[enr_full@result$Description %in% inter$a5]<-20 # all
enr_full@result$qvalue[enr_full@result$Description %in% inter$a5]<-20 # all
fname_net<-paste0(out_compare,'network_MOFA_single_',use_mofa_s, 
                  int_params ,N_EMAP,'.png')
create_multi_emap(enr_full, N_EMAP=100, fname_net=fname_net)
unique(inter$a7)




create_multi_emap<-function(enr_full, N_EMAP=200,fname_net, title){
  x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)
  as.character(enr_full@result$Description)
  
  
  
  emapplot(x2, showCategory=N_EMAP,
           cex_category=1, 
           cex_label_category=0.5, 
           color='qvalue', title=title,
  )
  
  ggsave(fname_net, width=15, height=15, dpi=300)
}



inter<-inter_mofa_union
unique(inter$a1)
#### Mofa vS single combination 
enr_full@result$Description<-enr_full@result$ID
enr_full@result$qvalue[enr_full@result$Description %in% inter$a3]<-10 # mirnas only

'cell death' %in% unique_mofa


create_multi_emap(enr_full, title='Mofa vs Union')
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


