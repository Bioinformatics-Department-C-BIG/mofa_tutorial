

###########
# create multi omics emap plot 
library('multienrichjam')
library('clusterProfiler')
library('enrichplot')

colnames(merged_factors_mofa)
colnames(enrich_rna)

merged_results<-merged_factors_mofa
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
merged_results$geneColname<-enrich_rna[match(merged_results$Description, enrich_rna$Description),]$core_enrichment
enr_full <- multienrichjam::enrichDF2enrichResult(as.data.frame(merged_results),
                                                  keyColname =  'Description',
                                                  geneColname ='geneColname',
                                                  pvalueColname = 'fish',
                                                  descriptionColname = 'Description',
                                                  qvalue='',
                                                  pvalueCutoff = 0.05)

enr_full@result$pvalue




# run_enrichment_plots(enr_full,results_file = 'test.txt')
enr_full@result$ID
enr_full@result$Description<-enr_full@result$ID
enr_full@result$qvalue[enr_full@result$Description %in% inter$a2]<-10 # mirnas only
enr_full@result$qvalue[enr_full@result$Description %in% inter$a10]<-5 # rna only
enr_full@result$qvalue[enr_full@result$Description %in% inter$a3]<-1 # mofa unique
enr_full@result$qvalue[enr_full@result$Description %in% inter$a5]<-20 # all


enr_full@result$qvalue


x2 <- pairwise_termsim(enr_full, showCategory = 180)
as.character(enr_full@result$Description)



emapplot(x2, showCategory=180,
         cex_category=1, 
         cex_label_category=0.5, 
         color='qvalue')


#### COMPARE SINGLE COMBINATION TO MOFA COMBINATION ####

cor_t<-0.15
outdir
#install.packages('arules')
library(arules)
library(igraph)
## create an emap as an igrpah object
ig_plot<-enrichMapJam(x2, n=120)
plot(ig_plot)



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


