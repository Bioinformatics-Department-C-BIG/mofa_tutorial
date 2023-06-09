

#script_dir<-dirname(rstudioapi::getSourceEditorContext()$path)
source(paste0(script_dir,'ppmi/setup_os.R'))


#### Metascripts 
#('UpSetR')
library('UpSetR')
library('dplyr')

library('VennDiagram')
library(grid)
#source(paste0(script_dir, 'ppmi/utils.R'))
#source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))

process_mirnas=FALSE

### Table of samples from all visits 
#### For each modality separately
out_compare<-'ppmi/plots/single/compare/'


cohort_cors
concatenate_pvals<-function(enrich_proteins,enrich_rna, enrich_mirnas=FALSE, pval_to_use='p.adjust') {
  #' supply the three enrichment datasets 
  #' chooose pval
  #' 
  #' 
  enrich_proteins_pvals<-enrich_proteins[, c(pval_to_use, 'Description')]
  enrich_rna_pvals<-enrich_rna[, c(pval_to_use, 'Description')]
  
  hist(enrich_proteins[,pval_to_use ])
  hist(enrich_rna[,pval_to_use ])
  # all=TRUE does not work well many mirna high are coming up 
  merged_paths<-merge(enrich_proteins_pvals, enrich_rna_pvals, by='Description', all=TRUE);dim(merged_paths)
  merged_paths<-merge(enrich_proteins_pvals, enrich_rna_pvals, by='Description');dim(merged_paths)
  if (add_mirs){
    enrich_mirna_pvals<-enrich_mirnas[, c(pval_to_use, 'Description')]
    hist(enrich_mirnas[,pval_to_use ])
    
    ## PROBLEM here if add all then something that is only in mirs it will come up 
    merged_paths<-merge(merged_paths,enrich_mirna_pvals, by='Description')
  }
  return(merged_paths)
  
}


get_combined_pvalue=function(merged_paths, pmethod='stouffer', weights=c(1,1,0.5)){
  library(metapod)
  
  p1<-merged_paths[,2]; length(p1)
  p2<-merged_paths[,3];length(p2)
  
  
  
  hist(p1)
  hist(p2)
  if (add_mirs){
    p3<-merged_paths[,4];length(p3)
    fish <- metapod::combineParallelPValues(list(p1, p2, p3),
                                            method=pmethod, 
                                            weights = weights)$p.value
    
  }else{
    fish <- metapod::combineParallelPValues(list(p1, p2),method=pmethod)$p.value
    
  }
  
  ### Add the combined to the original frame 
  merged_paths_fish<-cbind(merged_paths, fish)
  # and order
  merged_paths_fish<-merged_paths_fish[order(merged_paths_fish$fish),]
  
  ### Write significant results 
  
  cols_ch=colnames(merged_paths_fish)
  cols_ch=cols_ch[cols_ch!='Description']
  
  if (pval_to_use=='pvalue'){
    # adjust afterwards 
    merged_paths_fish_adj=merged_paths_fish
    merged_paths_fish_adj[cols_ch]=merged_paths_fish_adj[cols_ch]*dim(merged_paths_fish_adj)[1]
    merged_paths_fish_adj[cols_ch]=replace(merged_paths_fish_adj[cols_ch], merged_paths_fish_adj[cols_ch]>1, 1)
    merged_paths_fish=merged_paths_fish_adj
  }

  
  #if (add_mirs){
  #  merged_paths_fish_log$p.adjust<--log10(merged_paths_fish$p.adjust)
  #  
  #}
  if (use_mofa){
    write.csv(merged_paths_fish,paste0( merged_path_file2, '.csv'))
    
  }
  write.csv(merged_paths_fish,paste0( merged_path_file, '.csv'))
  
  return(merged_paths_fish)
  
}


create_combination_plot<-function(use_mofa, merged_paths_fish, combined_p_thresh=0.05, text_add=''){
  #### 
  #' @param  merged_paths_fish_log_sig
  #' @param combined_p_thresh default pvalue threshold=0.05
  #' @text_add add text to description
  #' 
  #' 
  #' 
  #combined_p_thresh=0.05
  # text_add=''
  # test
  #merged_paths_fish= mofa_only_all
  cols_ch=colnames(merged_paths_fish)[-1]
  merged_paths_fish_log=merged_paths_fish
  merged_paths_fish_log[cols_ch]=lapply(merged_paths_fish[cols_ch], 
                                        function(x){-log10(x)})
  
  
  merged_paths_fish_log_sig<-merged_paths_fish_log[merged_paths_fish_log$fish>-log10(combined_p_thresh),]
  
  merged_paths_fish_sig_melt<-reshape2::melt(merged_paths_fish_log_sig[1:Npaths,])
  merged_paths_fish_sig_melt
  if (add_mirs){
    labels<-c('proteomics', 'RNA', 'miRNA',  'Stouffer\'s')
  }else{ 
    labels<-c('proteomics', 'RNA', 'Stouffer\'s')}
  
  
  merged_paths_fish_sig_melt_fish<-merged_paths_fish_sig_melt[merged_paths_fish_sig_melt$variable=='fish',]
  plot_all=TRUE
  merged_to_plot<-merged_paths_fish_sig_melt
  
  if (plot_all){
    labels_p<-labels; width_p<-0.7
  }else{
    labels_p=labels[4];width_p<-0.7/4
  }
  
  NN=25
  nsig<-dim(merged_paths_fish_sig)[1]
  Npaths=ifelse(nsig<NN,nsig,NN)
  
  ## print the total number of paths 
  un_paths_0.05<-length(unique(merged_paths_fish_log$Description[merged_paths_fish_log$fish>-log10(0.05)]))
  un_paths_0.01<-length(unique(merged_paths_fish_log$Description[merged_paths_fish_log$fish>-log10(0.01)]))
  text_p<-paste0('\n p-adj.< ', 0.05,': ', un_paths_0.05, ' pathways', 
                 '\n p-adj.< ', 0.01,': ', un_paths_0.01, ' pathways \n' )
  
  if (use_mofa){
    
    text_p<-paste0(text_p,
                   'weighted: ', run_weighted, ', pval: ', pval_to_use, 
                   '\n mofa settings: ', TOP_PN,' ', TOP_GN,' ',TOP_MN)
    
    title_p<-paste( 'Factor ', sel_factors[fn],
                    ', cohort cor: ',  round(cohort_cors[sel_factors[fn]], digits=2), 
                    # '\n prot, RNA, miRNA \n', paste(round(weights_Var, digits=2),sep=' ',collapse=', '), 
                    '\n prot, RNA, miRNA \n', paste(round(vars_by_mod[c(3,2,1)], digits=2),sep=' ',collapse=', '))
    
  }else{
    title_p=''
  }
  
  title_p=paste0(title_p,'\n Weights: ', paste(round(weights_Var, digits=1), collapse=', '))
  mir_enrich_p_all<-ggplot(merged_to_plot,
                           aes( x=reorder(Description, value),
                                y=value, fill=variable))+
    #geom_bar(position='dodge', stat='identity', width=0.7)+
    geom_bar(position='dodge', stat='identity', width=width_p)+
    
    theme(axis.title.y=element_blank(), 
          axis.text.y= element_text(size=15,color='black' ), 
          axis.text.x= element_text(size=15,color='black' ),
          legend.title =element_blank(), 
          legend.position="top",
          plot.title =element_text(size=15, color='black'), 
          legend.text = element_text(size=10), 
          plot.caption= element_text(hjust = 0, face = "italic", size=15))+
    geom_hline(yintercept = -log10(combined_p_thresh), linetype='dashed')+
    labs(y='-log10pvalue', caption=text_p, 
         title =title_p)+
    scale_fill_discrete(labels=labels_p)+
    
    coord_flip()
  mir_enrich_p_all
  
  plot_params<-paste0(add_mirs,Npaths ,'_barplot_all',plot_all, combined_p_thresh,text_add)
  
  ggsave(paste0(merged_path_file, plot_params ,'.jpeg'),mir_enrich_p_all, dpi=300,
         width=floor(log10(Npaths))*3+7,height=log10(Npaths)*3+3 )
  if (use_mofa){
    ggsave(paste0(merged_path_file2, plot_params, '.jpeg'),mir_enrich_p_all, dpi=300,
           width=log10(Npaths)*3+7,height=log10(Npaths)*3+3)
  }
  mir_enrich_p_all
  
}



pval_to_use<-'pvalue'
pval_to_use<-'p.adjust'
add_mirs=TRUE
use_mofa=TRUE
v2=FALSE
if (v2){
  add_mirs=FALSE
  pval_to_use='p.adjust'
}
run_ORA=FALSE
pmethod<-'stouffer'

use_mofa
if (use_mofa){
  ### use all not just the significant p-value
  fn=4
  gse_mofa_rna=list1[[sel_factors[fn]]]
  cohort_cors[sel_factors[fn]]
  cohort_cors
  gse_mofa_prot=list_proteins[[sel_factors[fn]]]
  ## WARNING: RUN THE WHOLE enrich_mofa script to convert mirs 
  gse_mofa_mirs = list_mirs_enrich[[sel_factors[fn]]]
  enrich_rna<-gse_mofa_rna@result
  enrich_proteins=gse_mofa_prot@result
  enrich_mirnas=gse_mofa_mirs@result
  vars_by_mod<-vars_by_factor_all$r2_per_factor$group1[sel_factors[fn],]
  
  
  #### also weight by modality variance 
  t_var<-sum(vars_by_mod)
  mir_div=10
  weights_Var<-c(vars_by_mod['proteomics'], vars_by_mod['RNA'],vars_by_mod['miRNA']/mir_div )
  weights_Var<-weights_Var/sum(weights_Var)*100
  run_weighted=TRUE
  # also write to extra file
  use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
  
  merged_path_file2<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',  pvalueCutoff, pval_to_use,'_', run_ORA, pmethod,
                              'mofa_',  use_mofa_s , 'w_', run_weighted  )

  
}else{
  weights_Var=c(41,301,298)
  weights_Var=weights_Var/sum(weights_Var)*100
  run_weighted=TRUE
  enrich_rna<-enrich_rna_single
  enrich_proteins=enrich_proteins_single
  enrich_mirnas=enrich_mirnas_single
  
  
  use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
  
  
}


#### SET PARAMETERS

# TODO: ADD WEIGTS? IN FILE
merged_path_file<-paste0(out_compare, VISIT, '_',TISSUE,'_',  pvalueCutoff, pval_to_use,'_', run_ORA, pmethod,
                         'mofa_',  use_mofa_s  )



get_paths<-function(mode){
  prot_paths1<-read.csv(paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_',
                               T,  mode, '_enrichment_positive_pvals_no_f.csv' ))
  prot_paths2<-read.csv(paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_',
                               T,  mode, '_enrichment_positive_pvals_no_f.csv' ))
  prot_paths=rbind(prot_paths1,prot_paths2)
  colnames(prot_paths)<-c('X', 'Description', 'Factor', 'p.adjust')
  
  return(prot_paths)
}


if (use_mofa & v2){
  prot_paths<-get_paths('proteomics')
  rna_paths<-get_paths('RNA')
  prot_paths
  ### METHOD 2: PCGSE 
  fn=1
  ### use all not just the significant p-value
  cohort_cors[sel_factors[fn]]
  cohort_cors
  enrich_proteins=prot_paths[prot_paths$Factor==paste0('Factor', sel_factors[fn]),]
  enrich_rna<-rna_paths[rna_paths$Factor==paste0('Factor', sel_factors[fn]),]
  colnames(enrich_proteins)

  vars_by_mod<-vars_by_factor_all$r2_per_factor$group1[sel_factors[fn],]
  vars_by_mod
  vars_by_mod/sum(vars_by_mod) *100
  
}



################# ACTUALLY RUN THE combination #### 

if (add_mirs){
  merged_paths=concatenate_pvals(enrich_proteins,enrich_rna,enrich_mirnas, pval_to_use )
}else{
  merged_paths=concatenate_pvals(enrich_proteins=enrich_proteins,enrich_rna=enrich_rna, pval_to_use )
  
}


merged_paths_fish_res=get_combined_pvalue(merged_paths = merged_paths,weights= weights_Var)
merged_paths_fish=merged_paths_fish_res;dim(merged_paths_fish)
############################################
#### Create plots of combined p-values ####



create_combination_plot(use_mofa=use_mofa, merged_paths_fish=merged_paths_fish)
## question: what was NOT there before and is now? 




merged_path_file2
####### NOW MERGE THE RESULTS ####
f_pvals<-list()


###  LOAD SINGLE-COMBINATION ###


## what to actually input 
use_mofa=TRUE

### merge mofa 
for (fn in c(1,2,3,4)){
    use_mofa=TRUE;run_weighted=TRUE
    
    use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
    merged_path_file_load<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',  pvalueCutoff, pval_to_use,'_', run_ORA, pmethod,
                                  'mofa_',  use_mofa_s , 'w_', run_weighted  )
    
    f_pvals[[fn]]<-read.csv(paste0( merged_path_file_load, '.csv'), row.names = 1)
    
    
  }
  merged_factors_mofa<-do.call(rbind,f_pvals)
  merged_factors_mofa
  #### Merge factors together??? #### 
  dim(merged_factors_mofa[merged_factors_mofa$fish<0.05,])
  dim(merged_factors_mofa[merged_factors_mofa$fish<0.01,])

  
  
  ## LOAD SINGLE 
use_mofa=FALSE 
use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
merged_path_file_single<-paste0(out_compare, VISIT, '_',TISSUE,'_',  pvalueCutoff, pval_to_use,'_', run_ORA, pmethod,
                                  'mofa_',  use_mofa_s, 'w_', run_weighted  )
  
single_ps<-read.csv(paste0( merged_path_file_single, '.csv'), row.names = 1)
  




#### COMPARE TO SINGLE-UNION ####
#### option to compare MOFA or combination
### this is the single union of paths 

listInput_all_mods_single
combined_ps=single_ps

single_combination<-single_ps$Description[single_ps$fish< pvalueCutoff_sig]; 

listInput_single_combination=append(listInput_all_mods_single,
                                    list(single_combination), 
                                    )

inter<-calculate.overlap(listInput_single_combination)
write.csv(single_ps[combined_ps$Description %in% inter$a3,],
          paste0(out_compare, 'unique_single_comb_union_' ,use_mofa ,'.csv'))


names(listInput_single_combination)[4]<-'combination'
fname_venn=paste0(out_compare,'Single_union_vs_comb_all_modalities_', 
                  int_params ,'venn_mofa_',use_mofa, '.jpeg')

### def create a venn? 
venn_list=listInput_single_combination


create_venn<-function(venn_list, fname_venn){
  
  #######
  #' @param 
  #'
  #'
  myCol2 <- brewer.pal(length(venn_list), "Pastel2")
  venn.diagram(venn_list,
               # Circles
               lwd = 2, lty = 'blank', fill = myCol2, cex=2.5,cat.cex=1,
               filename =fname_venn, 
             
               output=FALSE)
}

create_venn(listInput_single_combination, fname_venn)

### create another list with all items 
merged_factors_mofa_ps<-merged_factors_mofa$Description[merged_factors_mofa$fish< pvalueCutoff_sig]; 

listInput_single_combination=append(listInput_all_mods_single,
                                    list(single_combination), 
)
# add also mofa
listInput_single_combination_all=append(listInput_single_combination,
                                        list(merged_factors_mofa_ps) 
)



inter<-calculate.overlap(listInput_single_combination_all)
write.csv(single_ps[combined_ps$Description %in% inter$a3,],
          paste0(out_compare, 'unique_single_comb_union_' ,use_mofa ,'.csv'))


names(listInput_single_combination_all)[4]<-'combination'
names(listInput_single_combination_all)[5]<-'mofa'

fname_venn=paste0(out_compare,'ALL_3', 
                  int_params ,'.jpeg')

### def create a venn? 
venn_list=listInput_single_combination_all
create_venn(listInput_single_combination_all, fname_venn)



### def create a venn? 






library(venneuler)
#install.packages('venneuler')

Lists <- listInput_single_combination_all  #put the word vectors into a list to supply lapply
items <- sort(unique(unlist(Lists)))   #put in alphabetical order
MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=length(Lists))  #make a matrix of 0's
colnames(MAT) <- names(Lists)
rownames(MAT) <- items
lapply(seq_along(Lists), function(i) {   #fill the matrix
  MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
})
dev.off()
MAT   #look at the results
v <- venneuler(MAT, quantities=TRUE)
plot(v, main=paste0('weights: ',paste0(round(weights_Var, digits=0), collapse = ', ')))







###########


#### COMPARE SINGLE COMBINATION TO MOFA COMBINATION ####

cor_t<-0.15
outdir

mofa_enrich_dir=paste0(outdir,'/enrichment/', mofa_params, TISSUE, 'ranked_list', cor_t)
mofa_enrich_file<-paste0(mofa_enrich_dir, '.csv')
all_ord_R<-read.csv(mofa_enrich_file, header=1)


### compare to merged 
p_final=0.05
all_ord_R<-merged_factors_mofa[merged_factors_mofa$fish<p_final,]
single_paths<-single_ps[single_ps$fish<p_final,]
dim(all_ord_R)

single_paths$Description<-gsub('-', ' ', tolower(single_paths$Description))
all_ord_R$Description<-gsub('-', ' ', tolower(all_ord_R$Description))

'cytokine mediated signaling pathway' %in% all_ord_R$Description
'cytokine mediated signaling pathway' %in% single_paths$Description

dim(single_paths)
dim(all_ord_R)
listInput_single_mofa<-list( mofa=unique(all_ord_R$Description), single=single_paths$Description)
mofa_overlap<-calculate.overlap(listInput_single_mofa)
unique_single<-mofa_overlap$a2[!(mofa_overlap$a2 %in%mofa_overlap$a3)]
unique_mofa<-mofa_overlap$a1[!(mofa_overlap$a1 %in%mofa_overlap$a3)]
mofa_overlap
unique_mofa
unique_single

mofa_overlap$a2
# common 
mofa_overlap$a3


myCol2 <- brewer.pal(4, "Pastel2")[1:length(listInput_single_mofa)]
venn.diagram(listInput_single_mofa,  # Circles
             lwd = 2, lty = 'blank',           fill = myCol2,
             cex=2.5, cat.cex=1,
             category.names=c( 'mofa\nintegration','single\ncombination'),
             filename = paste0(out_compare,'all_modalities_', 
                               int_params ,'venn_diagramm_mofa', cor_t, '.png'), 
             output=TRUE)



###### PLOT UNIQUE ONLY #### 
head(mofa_only_all)
mofa_only_all<-all_ord_R[all_ord_R$Description%in% unique_mofa,]
single_only_all<-single_ps[single_ps$Description%in% unique_single,]
mofa_only_all$Description
create_combination_plot(use_mofa=use_mofa, merged_paths_fish=mofa_only_all, combined_p_thresh=0.01, text_add='mofa_uniq')
create_combination_plot(use_mofa=use_mofa, merged_paths_fish=single_only_all, combined_p_thresh=0.01, text_add='single_uniq')


dim(mofa_only_all)
length(unique_mofa)

