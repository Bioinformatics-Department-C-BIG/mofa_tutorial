

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
enrich_mirnas


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
create_multi_emap(enr_full, N_EMAP=320)
unique(inter$a7)



inter<-inter_mofa_union
unique(inter$a1)
#### Mofa vS single combination 
enr_full@result$Description<-enr_full@result$ID
enr_full@result$qvalue[enr_full@result$Description %in% inter$a3]<-10 # mirnas only

'cell death' %in% unique_mofa


create_multi_emap(enr_full)
enr_full@result$qvalue


create_multi_emap<-function(enr_full, N_EMAP=200){
  x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)
  as.character(enr_full@result$Description)
  
  fname_net<-paste0(out_compare,'network_',use_mofa_s, 
                    int_params ,N_EMAP,'.png')
  
  emapplot(x2, showCategory=N_EMAP,
           cex_category=1, 
           cex_label_category=0.5, 
           color='qvalue', title('g')
            )
  
  ggsave(fname_net, width=15, height=15, dpi=300)
}

#### COMPARE SINGLE COMBINATION TO MOFA COMBINATION ####

cor_t<-0.15
outdir
#install.packages('arules')
library(arules)
library(igraph)
## create an emap as an igrpah object
N_EMAP=100
x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)

ig_plot<-enrichMapJam(x2, n=N_EMAP, 
                      )

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


