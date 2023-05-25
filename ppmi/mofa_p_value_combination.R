

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
  write.csv(merged_paths_fish,paste0( merged_path_file, '.csv'))
  
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
  merged_paths_fish_log=merged_paths_fish
  
  
  merged_paths_fish_log[cols_ch]=lapply(merged_paths_fish[cols_ch], 
                                        function(x){-log10(x)})
  
  #if (add_mirs){
  #  merged_paths_fish_log$p.adjust<--log10(merged_paths_fish$p.adjust)
  #  
  #}
  return(list(merged_paths_fish, merged_paths_fish_log))
  
}





#### Here is the input for the combination  ####
enrich_rna<-read.csv(enrich_rnas_file)
enrich_rna
dim(enrich_rna); dim(enrich_mirnas)
#enrich_mirnas<-read.csv(enrich_mirnas_file)
enrich_proteins<-read.csv(enrich_proteins_file)
enrich_mirnas<-read.csv(enrich_mirnas_file)
pval_to_use<-'p.adjust'
pval_to_use<-'pvalue'

add_mirs=TRUE
use_mofa=TRUE

if (use_mofa){
  
  ### use all not just the significant p-value
  fn=1
  gse_mofa_rna=list1[[sel_factors[fn]]]
  dim(gse_mofa_rna@result)
  
  gse_mofa_prot=list_proteins[[sel_factors[fn]]]
  ## WARNING: RUN THE WHOLE enrich_mofa script to convert mirs 
  gse_mofa_mirs = list_mirs_enrich[[sel_factors[fn]]]
  
  enrich_rna<-gse_mofa_rna@result
  enrich_proteins=gse_mofa_prot@result
  enrich_mirnas=gse_mofa_mirs@result
  
  vars_by_mod<-vars_by_factor_all$r2_per_factor$group1[sel_factors[fn],]
  vars_by_mod
}

vars_by_mod/sum(vars_by_mod) *100

################# ACTUALLY RUN THE combination #### 

### WHETHER TO ADD MIRNA INFO

if (add_mirs){
  merged_paths=concatenate_pvals(enrich_proteins,enrich_rna,enrich_mirnas, pval_to_use )
  
}else{
  merged_paths=concatenate_pvals(enrich_proteins,enrich_rna, pval_to_use )
  
}
sapply(list(enrich_rna,
            enrich_mirnas, 
            enrich_proteins), function(df){  
              length(which(df[, pval_to_use]<0.05))})

use_mofa
merged_paths

### Define files and parameters
use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
#merge(enrich_proteins)
pmethod<-'stouffer'
merged_path_file<-paste0(out_compare, VISIT, '_',TISSUE,'_',  pvalueCutoff, pval_to_use,'_', run_ORA, pmethod,
                         'mofa_',  use_mofa_s  )
merged_path_file2<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',  pvalueCutoff, pval_to_use,'_', run_ORA, pmethod,
                          'mofa_',  use_mofa_s  )


### Obtain the combined value and the log pvals
#  TODO: CAREFUL what happens if i dont give mirna
merged_paths_fish_res=get_combined_pvalue(merged_paths = merged_paths)
t_var<-sum(vars_by_mod)
weights_Var<-c(vars_by_mod['proteomics'], vars_by_mod['RNA'],vars_by_mod['miRNA']/30 )
#weights_Var=log(weights_Var*100)
weights_Var
merged_paths[merged_paths$Description=='inflammatory response' ,]
merged_paths_fish_res=get_combined_pvalue(merged_paths = merged_paths,weights= weights_Var)

#merged_paths_fish_res=get_combined_pvalue(merged_paths = merged_paths)
merged_paths_fish=merged_paths_fish_res[[1]];dim(merged_paths_fish)

merged_paths_fish_padj<-is.numeric(merged_paths_fish)/dim(merged_paths_fish)[1]
merged_paths_fish_log=merged_paths_fish_res[[2]];dim(merged_paths_fish_log)
head(merged_paths_fish)
merged_paths_fish[merged_paths_fish$Description=='inflammatory response' ,]

dim(merged_paths_fish)[1]
combined_p_thresh<-0.05
length(which(merged_paths$p.adjust.x<combined_p_thresh))
length(which(merged_paths$p.adjust.y<combined_p_thresh))
length(which(merged_paths$p.adjust<combined_p_thresh))

merged_paths_fish_sig<-merged_paths_fish[merged_paths_fish$fish<combined_p_thresh,]
merged_paths_fish_log_sig<-merged_paths_fish_log[merged_paths_fish_log$fish>-log10(combined_p_thresh),]
dim(merged_paths_fish_log_sig)

############################################
#### Create plots of combined p-values ####

Npaths=25

merged_paths_fish_sig_melt<-reshape2::melt(merged_paths_fish_log_sig[1:Npaths,], 
                                           varnames=c('modality','nlog10pvalue' ))

if (add_mirs){
  labels<-c('proteomics', 'RNA', 'miRNA',  'Stouffer\'s')
}else{ 
  labels<-c('proteomics', 'RNA', 'Stouffer\'s')}


merged_paths_fish_sig_melt_fish<-merged_paths_fish_sig_melt[merged_paths_fish_sig_melt$variable=='fish',]
merged_paths_fish_sig_melt_fish
plot_all=TRUE

if (plot_all){
  merged_to_plot<-merged_paths_fish_sig_melt
  labels_p<-labels; width_p<-0.7
}else{
  merged_to_plot<-merged_paths_fish_sig_melt_fish
  labels_p=labels[4];width_p<-0.7/4
}


merged_paths_fish_sig[merged_paths_fish_sig$Description=='inflammatory response',]
## print the total number of paths 
un_paths<-dim(merged_paths_fish_sig)[1]
text_p<-paste0('\n p-adj.< ', combined_p_thresh,': ', un_paths, ' pathways')

mir_enrich_p_all<-ggplot(merged_to_plot,
                         aes( x=reorder(Description, value),
                              y=value, fill=variable))+
  #geom_bar(position='dodge', stat='identity', width=0.7)+
  geom_bar(position='dodge', stat='identity', width=width_p)+
  
  theme(axis.title.y=element_blank(), 
        axis.text.y= element_text(size=15,color='black' ), 
        axis.text.x= element_text(size=15,color='black' ),
        legend.title =element_blank(), 
        legend.text = element_text(size=15), 
        plot.caption= element_text(hjust = 0, face = "italic", size=20))+
  labs(y='-log10pvalue', caption=text_p)+
  scale_fill_discrete(labels=labels_p)+
  
  coord_flip()
mir_enrich_p_all
ggsave(paste0(merged_path_file, add_mirs,Npaths ,'_barplot_all',plot_all, combined_p_thresh, '.jpeg'),mir_enrich_p_all, dpi=300,
       width=Npaths*0.55,height=Npaths/2.5 )

ggsave(paste0(merged_path_file2, add_mirs,Npaths ,'_barplot_all',plot_all, combined_p_thresh, '.jpeg'),mir_enrich_p_all, dpi=300,
       width=Npaths*0.55,height=Npaths/2.5 )


Npaths

