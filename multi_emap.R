
source(paste0(script_dir,'ppmi/mofa_p_value_combination.R'))
source(paste0(script_dir,'ppmi/utils.R'))

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

f_pvals[2]

merged_results<-f_pvals_merged



### TODO: where to obtain the gene names from??? 
# Now they are taken from RNA 
merged_results$geneColname<-enrich_rna[match(merged_results$Description, enrich_rna$Description),]$core_enrichment

enrich_res_df<-merged_results



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
create_multi_emap<-function(enr_full, N_EMAP=25,fname_net, title, color_by='qvalue', min=0.2){
          #'
          #'
          #'
          #'
          
          #N_EMAP=50
          #min=0.3
  
          options(ggrepel.max.overlaps = Inf)
          
          x2 <- pairwise_termsim(enr_full, showCategory = N_EMAP)
          
          
          
          
          
          
          color_by='pvalue'
          #jpeg(fname_net)
          em<-emapplot(x2, showCategory=N_EMAP,
                       cex_category=1, 
                       cex_label_category=0.5, 
                       color=color_by,
                       min =min, 
                          
          )
          em
          fn
          
          
  
}




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
        ## obtain the enrich RNA of the specific factor!!! 
        weights_Var=get_mofa_vars(factor=sel_factors[fn], adj_weights)
        
        
        factor1_paths$geneColname<-enrich_rna[match(factor1_paths$Description, enrich_rna$Description),]$core_enrichment
        which.max(weights_Var)
        if (which.max(weights_Var)=='miRNA'){
          factor1_paths$geneColname<-enrich_mirnas[match(factor1_paths$Description, enrich_mirnas$Description),]$geneID
          
        }else if(which.max(weights_Var)=='proteomics'){
          factor1_paths$geneColname<-enrich_proteins[match(factor1_paths$Description, enrich_proteins$Description),]$geneID
          
          
        }
        
        enrich_res_df<-factor1_paths

        enr_full<-create_enrich_result(enrich_res_df)
        
        
        # run_enrichment_plots(enr_full,results_file = 'test.txt')
        enr_full@result$Description<-enr_full@result$ID
        
        
        min_overlap_score=0.3
        #### EMAP by factor ####
        create_multi_emap(enr_full, N_EMAP=30, fname_net=paste0(merged_path_file_mofa, '_network', '.jpeg'), color_by = 'pvalue', 
                          min = min_overlap_score)
        


}

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
  
  
  ### ALL FACTOR PATHS 
  merged_results_sig=merged_results[merged_results$fish<0.05,]
  merged_results_sig=merged_results_sig[c('Description', 'fish')]
  ## obtain the enrich RNA of the specific factor!!! 

  ### Connect by gene OVERLAP 
  merged_results_sig$geneColname<-enrich_rna[match(merged_results_sig$Description, enrich_rna$Description),]$core_enrichment
  
  
  enrich_res_df<-merged_results_sig

  enr_full<-create_enrich_result(enrich_res_df)
  
  
  # run_enrichment_plots(enr_full,results_file = 'test.txt')
  enr_full@result$Description<-enr_full@result$ID
  
  
  min_overlap_score=0.2
  #### EMAP by factor ####
  create_multi_emap(enr_full, N_EMAP=400, fname_net=paste0(merged_path_file_mofa, '_network_ALLfs', '.jpeg'), color_by = 'pvalue', 
                    min_edge = min_overlap_score)
  
  
  















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

  
  
  ggsave(fname_net, width=9, height=9, 
         dpi = 300)  #ggsave(filename=fname_net, plot=last_plot(),width=15, height=15, dpi=300)
}



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



