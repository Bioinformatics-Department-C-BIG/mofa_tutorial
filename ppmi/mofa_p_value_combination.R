script_dir
#setwd(script_dir)

source(paste0('ppmi/setup_os.R'))

library('UpSetR')
library('dplyr')
library('ggplot2')
library('VennDiagram')
library(grid)
script_dir
source(paste0(script_dir, 'ppmi/utils.R'))
#source(paste0(script_dir,'ppmi/deseq_analysis_setup.R'))
## load single combination enrichment
#source(paste0(script_dir,'ppmi/load_single_combination_pathways.R'))
### load mofa enrichment 
source(paste0(script_dir,'ppmi/mofa_enrich.R'))


standardize_go_names<-function(descriptions){
  #'
  #'
  #'
  descriptions=gsub('-', ' ', tolower(descriptions))
  descriptions=gsub('_', ' ', tolower(descriptions))
  
  #descriptions=gsub('^[:alnum:]', '', tolower(descriptions))
  

  descriptions=gsub("\\'", '', tolower(descriptions))
  descriptions=gsub("\\,", '', tolower(descriptions))
  
  descriptions=gsub('\\(', '', tolower(descriptions))
  descriptions=gsub('\\)', '', tolower(descriptions))
  descriptions=gsub('\\/', '', tolower(descriptions))
  
  
  
  
  return(descriptions)
  
}





process_mirnas=FALSE

### Table of samples from all visits 
#### For each modality separately
out_compare<-'ppmi/plots/single/compare/'
out_compare<-paste0('ppmi/plots/single/compare/', VISIT, '/')
dir.create(out_compare)


cohort_cors
concatenate_pvals <- function(enrich_proteins,enrich_rna, enrich_mirnas=FALSE, pval_to_use='p.adjust') {
  #' supply the three enrichment datasets 
  #' chooose pval
  #' 
  #' 
  enrich_proteins_pvals<-as.data.frame(enrich_proteins[c(pval_to_use, 'Description')])
  enrich_rna_pvals<-as.data.frame(enrich_rna[c(pval_to_use, 'Description')])
  
  #hist(enrich_proteins[,pval_to_use ])
  #hist(enrich_rna[,pval_to_use ])
  # all=TRUE does not work well many mirna high are coming up 
  merged_paths<-merge(enrich_proteins_pvals, enrich_rna_pvals, by='Description', all=TRUE);
  dim(merged_paths)
  #merged_paths<-merge(enrich_proteins_pvals, enrich_rna_pvals, by='Description');
  dim(merged_paths)
  ## TODO: double check the merging here -- maybe there are dashes and different formats 
  if (add_mirs){
    enrich_mirna_pvals<-enrich_mirnas[, c(pval_to_use, 'Description')]
    #hist(enrich_mirnas[,pval_to_use ])
    
    ## PROBLEM here if add all then something that is only in mirs it will come up 
    merged_paths<-merge(merged_paths,enrich_mirna_pvals, by='Description')
  }
  return(merged_paths)
  
}

#install.packages('sgof')
library(sgof)
get_combined_pvalue=function(merged_paths, pmethod='stouffer', weights_Var=c(1,1,0.5), merged_path_file, pval_to_use='pvalue'){
  #'
  #' Combine and save 
  
  #View(merged_paths_fish[grep('oxidative', merged_paths_fish$Description),])
  #merged_paths[grep('oxidative phosphorylation', merged_paths$Description),]$pvalue.x=0.999999
  #merged_paths[grep('oxidative phosphorylation', merged_paths$Description),]
  
  
  library(metapod)
  merged_paths[is.na(merged_paths)]<-0.999
  p1<-merged_paths[,2]; length(p1)# proteins 
  p2<-merged_paths[,3];length(p2)
  p3<-merged_paths[,3];length(p3)
  #weights=weights_Var
  print(paste0('Weights ', weights_Var))
  list_total=list()
  weights_Var_ge1=c()
  if (weights_Var['proteomics']>1){
    list_total<-append(list_total, list(p1))
    weights_Var_ge1=c(weights_Var['proteomics'])
  }  
  if (weights_Var['RNA']>1){
    list_total<-append(list_total, list(p2))
    weights_Var_ge1=c(weights_Var_ge1,weights_Var['RNA'])
    
  }
  if (weights_Var['miRNA']>1){
    list_total<-append(list_total, list(p3))
    weights_Var_ge1=c(weights_Var_ge1,weights_Var['miRNA'])
    
  }
  
  print(paste('weights ',round(weights_Var_ge1)))
  length(list_total)
  #weights_Var=c(1/500, 5/100, 5/100)
  hist(p1)
  hist(p2)
  if (add_mirs){
    #p3<-merged_paths[,4];length(p3)
    fish <- metapod::combineParallelPValues(list_total,
                                            method=pmethod, 
                                            weights = weights_Var_ge1)$p.value
    
  }else{
    fish <- metapod::combineParallelPValues(list(p1, p2),method=pmethod)$p.value
    
  }
  
  ### Add the combined to the original frame 
  merged_paths_fish<-cbind(merged_paths, fish)
  head(merged_paths_fish)
  # and order
  merged_paths_fish<-merged_paths_fish[order(merged_paths_fish$fish),]

    ### Write significant results 
  
  cols_ch=colnames(merged_paths_fish)
  cols_ch=cols_ch[cols_ch!='Description']
  
  if (pval_to_use=='pvalue'){
    # adjust afterwards 
    #merged_paths_fish_adj=merged_paths_fish
    merged_paths_fish$fish_pvalue<-merged_paths_fish$fish
    #print(paste('Adjusted p values using: ', dim(merged_paths_fish_adj)[1]))
    
  #  merged_paths_fish_adj[cols_ch]=merged_paths_fish_adj[cols_ch]*dim(merged_paths_fish_adj)[1]
   # merged_paths_fish_adj[cols_ch]=replace(merged_paths_fish_adj[cols_ch], merged_paths_fish_adj[cols_ch]>1, 1)
    #merged_paths_fish=merged_paths_fish_adj
    ## apply benjamini-holberg
    #length(which(merged_paths_fish$fish	<0.05))
    
    merged_paths_fish$fish<-BH(merged_paths_fish$fish, alpha = 0.05)$Adjusted.pvalues	

  }

  
  #if (add_mirs){
  #  merged_paths_fish_log$p.adjust<--log10(merged_paths_fish$p.adjust)
  #  
  #}

  write.csv(merged_paths_fish,paste0( merged_path_file, '.csv'))
  which(merged_paths_fish$fish<0.05)
  print(paste('Saving ',merged_path_file ))
  return(merged_paths_fish)
  
}


create_combination_plot<-function(use_mofa, merged_paths_fish, combined_p_thresh=0.05, text_add='', title_p='', merged_path_file, plot_all=TRUE){
  ####  Creates a plot from the combined result from Stouffer's 
  #' @param  merged_paths_fish_log_sig expects order of columns as:  proteins, rnas, mirnas 
  #' @param combined_p_thresh default pvalue threshold=0.05
  #' @text_add add text to description
  #' 
  #' 
  #' 
  #combined_p_thresh=0.05
  # text_add=''
  # test
  cols_ch=colnames(merged_paths_fish)[-1]
  merged_paths_fish_log=merged_paths_fish
  merged_paths_fish_log[cols_ch]=lapply(merged_paths_fish[cols_ch], 
                                        function(x){-log10(x)})
  
  
  merged_paths_fish_log_sig<-merged_paths_fish_log[merged_paths_fish_log$fish>-log10(combined_p_thresh),]
  
  #merged_paths_fish_log_sig=merged_paths_fish_log_sig[merged_paths_fish_log_sig$fish,]
  
  
  ### FILTER THE FIRST NPATHS HERE
  NN=25
  nsig<-dim(merged_paths_fish_log_sig)[1];nsig
  Npaths=ifelse(nsig<NN,nsig,NN)
  
  
  merged_paths_fish_sig_melt<-reshape2::melt(merged_paths_fish_log_sig[1:Npaths,])
  merged_paths_fish_sig_melt
  if (add_mirs){
    labels<-c('proteomics', 'RNA', 'miRNA',  'Stouffer\'s')
  }else{ 
    labels<-c('proteomics', 'RNA', 'Stouffer\'s')}
  
  # filter only fish
  merged_paths_fish_sig_melt_fish<-merged_paths_fish_sig_melt[merged_paths_fish_sig_melt$variable=='fish',]
  plot_all=TRUE
  merged_to_plot<-merged_paths_fish_sig_melt

  if (plot_all){
    labels_p<-labels; width_p<-0.7
    
  }else{
    labels_p=labels[4];width_p<-0.7/4
    
  }
  

  
  
  ## ORDER BY DECREASING FISH
  ord_paths<-merged_paths_fish_sig_melt_fish$Description[order(merged_paths_fish_sig_melt_fish$value)]
  merged_paths_fish_sig=merged_paths_fish[merged_paths_fish$fish<0.05,]
 
  
  ## print the total number of paths 
  un_paths_0.05<-length(unique(merged_paths_fish_log$Description[merged_paths_fish_log$fish>-log10(0.05)]))
  un_paths_0.01<-length(unique(merged_paths_fish_log$Description[merged_paths_fish_log$fish>-log10(0.01)]))
  text_p<-paste0('\n p-adj.< ', 0.05,': ', un_paths_0.05, ' pathways', 
                 '\n p-adj.< ', 0.01,': ', un_paths_0.01, ' pathways \n' )
  
  
  merged_to_plot$Description=factor(merged_to_plot$Description, levels=unique(ord_paths))
  levels(merged_to_plot$Description)
    

    
  mir_enrich_p_all<-ggplot(merged_to_plot,
                          # aes( x=reorder(Description, value),
                          aes( x=Description,
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
  

}


VISIT
pval_to_use<-'pvalue'
pval_to_use<-'p.adjust'
pval_to_use<-'pvalue'

add_mirs=TRUE
use_mofa=TRUE
v2=FALSE
if (v2){
  add_mirs=FALSE
  pval_to_use='p.adjust'
}
run_ORA=FALSE
pmethod<-'stouffer'

get_mofa_paths_and_weights<-function(factor){
        #'
        #'For each mofa factor extract the list of paths and the weights 
        #' @param fn factor value to obtain
        #' @return      
        
        gse_mofa_rna=list1[[factor]]
        gse_mofa_prot=list_proteins[[factor]]
        gse_mofa_mirs = list_mirs_enrich[[factor]]
        
        enrich_rna<-gse_mofa_rna@result
        enrich_proteins=gse_mofa_prot@result
        enrich_mirnas=gse_mofa_mirs@result
        
        
        return(list(enrich_rna, enrich_proteins, enrich_mirnas))
        
        
  }

        get_mofa_vars<-function(factor,adj_weights){
          #'
          #'
                vars_by_mod<-vars_by_factor_all$r2_per_factor$group1[factor,]
                
                
                #### also weight by modality variance 
                t_var<-sum(vars_by_mod)
                weights_Var1<-c(vars_by_mod['proteomics'], vars_by_mod['RNA'],vars_by_mod['miRNA'] )
                weights_Var1<-weights_Var1/sum(weights_Var1)*100
                weights_Var2<-adj_weights/sum(adj_weights)*100
                weights_Var=weights_Var1*weights_Var2
                weights_Var=weights_Var/sum(weights_Var)*100
                print(weights_Var)
                return(weights_Var)
          
        }
        
get_combination_settings<-function(weights_var,adj_weights=c(1,1,1), use_mofa=FALSE, fn=NULL, vars_by_mod=NULL){
      #'
      #' 
      #'
      #'      

        if (use_mofa){
          
          #  text_p<-paste0(text_p,
          #                ', pval: ', pval_to_use, 
          #               '\n mofa settings: ', TOP_PN,' ', TOP_GN,' ',TOP_MN)
          
          title_p<-paste( 'Factor ', sel_factors[fn],
                          ', cohort cor: ',  round(cohort_cors[sel_factors[fn]], digits=2), 
                          # '\n prot, RNA, miRNA \n', paste(round(weights_Var, digits=2),sep=' ',collapse=', '), 
                          '\n prot, RNA, miRNA \n weights adj.: ', paste(round(adj_weights, digits=2),sep=' ',collapse=', '))
          
        }else{
          title_p=''
        }
        
        title_p=paste0(title_p,'\n Weights: ', paste(round(weights_var, digits=1), collapse=', '))
        
        return(title_p)
  }      
        
        
### use all not just the significant p-value
        ## run for all mofa factors and concatenate? 
use_mofa=TRUE;run_weighted=TRUE
f_pvals<-list()
pval_to_use<-'p.adjust'
pval_to_use<-'pvalue'

v2=FALSE

## selected factor  indices in numbers 
fns=c(1:length(sel_factors))

        for (fn in fns){
                #fn=2
                print(sel_factors[fn])
                list_all<-get_mofa_paths_and_weights(factor=sel_factors[fn])
                # also write to extra file
                
                enrich_rna= list_all[[1]]
                
                if (v2){
                  res_merged_l<-read.csv(paste0(outdir,'/enrichment/' ,gsub('\\:', '_', subcategory), '_', T, mode, '_enrichment.csv' ))
                  res_merged_l_f<-res_merged_l[res_merged_l$Factor ==names(sel_factors[fn]),]
                }
                
                res_merged_l_f$Description<-standardize_go_names(res_merged_l_f$Description)
                res_merged_l_f$Description<-gsub('gobp ','',res_merged_l_f$Description)
                res_merged_l_f$pvalue<-res_merged_l_f$pvalue_min
                
                #View(res_merged_l_f[grep('oxidative', res_merged_l_f$Description),])
                
                enrich_proteins = list_all[[2]]
                enrich_mirnas = list_all[[3]]
                #View(enrich_mirnas[grep('oxidative', enrich_mirnas$Description),])
                
                adj_weights=c(41,301,298)
                adj_weights=c(1,1,1)
                
                weights_Var=get_mofa_vars(factor=sel_factors[fn], adj_weights)
                weights_Var
                print(paste(sel_factors[fn]))
                print(weights_Var, digits=2)
                use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
                merged_path_file_mofa<-paste0(outdir, '/enrichment/', VISIT, '_',TISSUE,'_',run_ORA, pmethod,
                                              'mofa_',  use_mofa_s   )
                
                title_p=get_combination_settings(weights_var=weights_Var, use_mofa=use_mofa,fn=fn, adj_weights=adj_weights)
                
               
                
                ################# ACTUALLY RUN THE combination #### 
                cors_pearson_l
                merged_paths=concatenate_pvals(enrich_proteins=enrich_proteins,
                                               enrich_rna=enrich_rna,enrich_mirnas=enrich_mirnas, pval_to_use=pval_to_use )
                
                if (v2){
                  enrich_proteins$Description=standardize_go_names(enrich_proteins$Description)
                  enrich_mirnas$Description=standardize_go_names(enrich_mirnas$Description)
                  
                  merged_paths=concatenate_pvals(enrich_proteins=enrich_proteins,
                                                 enrich_rna=res_merged_l_f,enrich_mirnas=enrich_mirnas, pval_to_use=pval_to_use )
                  merged_paths
                  }
                #merged_paths[merged_paths==1]=0.999999999999999
                #View(enrich_proteins[grep('oxidative', enrich_proteins$Description),])
                
                #View(merged_paths[grep('oxidative', merged_paths$Description),])
                
                merged_paths_fish_res=get_combined_pvalue(merged_paths = merged_paths,weights_Var=weights_Var, 
                                                          merged_path_file=merged_path_file_mofa, pval_to_use=pval_to_use)
                
                merged_paths_fish=merged_paths_fish_res;dim(merged_paths_fish)
                which(merged_paths_fish$fish<0.05)
                
                                #### Create plots of combined p-values ####
                
                # TODO: INVESTIGATE merged combination plot crashes at factor 3 
                #create_combination_plot(use_mofa=use_mofa, merged_paths_fish=merged_paths_fish,title_p=title_p, 
                #                        merged_path_file=merged_path_file_mofa )
                ## question: what was NOT there before and is now? 
                f_pvals[[fn]]<-merged_paths_fish
                
               # f_pvals[[fn]]<-read.csv(paste0( merged_path_file_mofa, '.csv'), row.names = 1)
               # lapply(f_pvals,dim)
                print(lapply(f_pvals, function(x){length(which(x$fish<0.05))}))
                
                
                
                
        }     


## how many significant for each factor       
lapply(f_pvals, function(x){length(which(x$fish<0.05))}) 
lapply(f_pvals, function(x){which(x$fish<0.05)}) 



#merged_factors_mofa[merged_factors_mofa$Description =='cell cycle',]

merged_factors_mofa<-do.call(rbind,f_pvals)
#### Merge factors together??? #### 
merged_factors_mofa_sig<-merged_factors_mofa[merged_factors_mofa$fish<0.05,]
merged_factors_mofa_sig$Description<-gsub('-', ' ', tolower(merged_factors_mofa_sig$Description))
unique(length(merged_factors_mofa_sig$Description))



#merged_factors_mofa_sig_uniq<-rank_mofa_paths(merged_factors_mofa_sig)
#head(merged_factors_mofa_sig_uniq)
## what to actually input 

### merge mofa 


use_mofa=FALSE
### single combination MODE ##### 
single_weights_Var=c(1,1,1)
single_weights_Var=c(41,301,298)

  #weights_Var=weights_Var/sum(weights_Var)*100
  run_weighted=TRUE
  enrich_rna<-enrich_rna_single
  enrich_proteins=enrich_proteins_single
  enrich_mirnas=enrich_mirnas_single

  use_mofa_s=ifelse(use_mofa, paste0('_',sel_factors[fn]),use_mofa )
  
  merged_path_file_single<-paste0(out_compare, VISIT, '_',TISSUE,'_', run_ORA, pmethod,
                           'mofa_',  use_mofa_s  )
  
  
  ################# ACTUALLY RUN THE combination #### 
  title_p=get_combination_settings(weights_var=single_weights_Var, use_mofa=FALSE,fn=fn, adj_weights=adj_weights)
  title_p
  merged_paths=concatenate_pvals(enrich_proteins=enrich_proteins,enrich_rna=enrich_rna,enrich_mirnas=enrich_mirnas, pval_to_use=pval_to_use )
  merged_paths_fish_res=get_combined_pvalue(merged_paths = merged_paths,weights= single_weights_Var,merged_path_file= merged_path_file_single, 
                                            pval_to_use=pval_to_use)
  merged_paths_fish=merged_paths_fish_res;dim(merged_paths_fish)
  #### Create plots of combined p-values ####
  create_combination_plot(use_mofa=use_mofa, merged_paths_fish=merged_paths_fish,title_p=title_p, merged_path_file=merged_path_file_single   )
  
### MAIN OUTPUT: 
  #$merged_path_file_single

#### SET PARAMETERS



  







####### NOW MERGE THE RESULTS ####

  
  
#### LOAD SINGLE COMBINATION ####
single_ps<-read.csv(paste0( merged_path_file_single, '.csv'), row.names = 1)
single_ps_sig<-single_ps$Description[single_ps$fish< pvalueCutoff_sig]; 
single_ps_sig<-gsub('-', ' ', tolower(single_ps_sig))
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


use_mofa=TRUE
library('RColorBrewer')
length(venn_list)
create_venn<-function(venn_list, fname_venn, main){
  
  #######
  #' @param 
  #'
  #'
  myCol2 <- brewer.pal(length(venn_list), "Pastel2")[1:length(venn_list)]
  venn.diagram(venn_list,
               # Circles
               lwd = 2, lty = 'blank', fill = myCol2, cex=2.5,cat.cex=1.5,
               filename =fname_venn, 
               main=main,
               
               output=FALSE)
}




library(venneuler)
#install.packages('venneuler')

#Lists <- listInput_single_combination_all  #put the word vectors into a list to supply lapply

create_venneuler<-function(Lists,fname){
  
  
  items <- sort(unique(unlist(Lists)))   #put in alphabetical order
  MAT <- matrix(rep(0, length(items)*length(Lists)), ncol=length(Lists))  #make a matrix of 0's
  colnames(MAT) <- names(Lists)
  rownames(MAT) <- items
  lapply(seq_along(Lists), function(i) {   #fill the matrix
    MAT[items %in% Lists[[i]], i] <<- table(Lists[[i]])
  })
  v <- venneuler(MAT, quantities=TRUE)
  
  jpeg(fname, res=300, height=1000, width=1000)
  plot(v, cex=1.1)
  dev.off()
  
}


#### COMPARE TO SINGLE MODALITIES SEPARATELY-UNION ####
#### option to compare MOFA or combination
### this is the single union of paths 

listInput_all_mods_single
use_mofa=TRUE
### CHOOSE WHICH OF THE TWO TO PLOT HERE 
if (use_mofa){
    combined_ps<-merged_factors_mofa_sig$Description
    title='MOFA vs Single'
  }else{
    combined_ps<-single_ps_sig
    title='Single combination vs Single'
  }
# combined_ps=merged_factors_mofa; dim(merged_factors_mofa)
listInput_all_mods_single$rna

listInput_single_combination=append(listInput_all_mods_single[c(1,3)],
                                    list(combined_ps))

inter_mofa_single<-calculate.overlap(listInput_single_combination)

names(listInput_single_combination)[length(listInput_single_combination)]<-'combination'
fname_venn=paste0(out_compare,'Single_union_vs_comb_all_modalities_', 
                  int_params ,'venn_mofa_',use_mofa, '.jpeg')

### def create a venn? 
venn_list=listInput_single_combination
create_venn(listInput_single_combination, fname_venn, main=title)

mofa_rna_common
fname_venn
##########################
#########################




#### Compare mofa to rna only! ####
listInput_all_mods_single$rna
merged_factors_mofa_sig
listInput_MOFA_RNAonly=append(list(listInput_all_mods_single$rna),
                                    list(unique(merged_factors_mofa_sig$Description)))
names(listInput_MOFA_RNAonly)<-c('rna', 'mofa')
fname_venn=paste0(out_compare,'MOFAvsRNA_', 
                  int_params ,'venn_mofa_',use_mofa, '.jpeg')

### def create a venn? 
venn_list=listInput_MOFA_RNAonly
mofa_rna_common<-intersect(listInput_MOFA_RNAonly$rna,listInput_MOFA_RNAonly$mofa )
create_venn(listInput_MOFA_RNAonly, fname_venn, main='MOFA vs RNA only')

mofa_unique<-listInput_MOFA_RNAonly$mofa[!listInput_MOFA_RNAonly$mofa %in%mofa_rna_common]
rna_unique<-listInput_MOFA_RNAonly$rna[!listInput_MOFA_RNAonly$rna %in%mofa_rna_common]


#### GET PARENTS  ### Exclude unique mofa parents and children from the really unique ones 
go_ids$TERM[grep('small gtpase', standardize_go_names(go_ids$TERM))]

mofa_unique_id<-go_ids$GOID[match(mofa_unique, standardize_go_names(go_ids$TERM))]
rna_unique_id<-go_ids$GOID[match(rna_unique, standardize_go_names(go_ids$TERM))]
go_ids_rna<-go_ids$GOID[grep(all_children[1], go_ids$TERM)]


## For each unique mofa path--> only set it as unique IF and ONLY if it does not have children in rna unique
all_children<-lapply(mofa_unique_id, function(x){get_child_nodes(x)})
all_parents<-lapply(mofa_unique_id, function(x){ps<-get_parent_nodes(x)
                                                ps<-ps[ps$distance<3,]
                                                return(ps)})


unique(rna_unique_id)

### get also parents of rna unique
#parents_rna_unique<-lapply(rna_unique_id[!is.na(rna_unique_id)], function(x){ps<-get_parent_nodes(x)
#                                                      ps<-ps[ps$distance<3,]
#                                                      return(ps)})
#children_rna_unique<-lapply(rna_unique_id[!is.na(rna_unique_id)], function(x){get_child_nodes(x)})

### Which of the children overlap with 
res_all<-unlist(lapply(all_children, function(df){
  res<-any(  df$child_go_id %in% c(rna_unique_id )    )
  return( res )
}))

res_all_ps<-unlist(lapply(all_parents, function(df){
  res<-any(  df$parent_go_id %in% c(rna_unique_id  ))
  return( res )
}))


sum((res_all))
mofa_unique_nic<-mofa_unique[!(res_all|res_all_ps) ]
mofa_unique_nic

mofa_unique_nic
listInput_MOFA_RNAonly$rna[grep( 'ras protein',standardize_go_names(listInput_MOFA_RNAonly$rna))]
listInput_MOFA_RNAonly$mofa[grep( 'ras protein',listInput_MOFA_RNAonly$mofa)]














############## Test mofa uniqueness ##############################



### create another list with all THREE items 
merged_factors_mofa_ps<-merged_factors_mofa_sig$Description
listInput_all_mods_single=lapply(listInput_all_mods_single, function(x) {gsub('-', ' ', tolower(x)) }) 
listInput_single_combination=append(listInput_all_mods_single[c(1,3)],
                                    list(single_ps_sig) 
)
# add also mofa
listInput_single_combination_all=append(listInput_single_combination,
                                        list(merged_factors_mofa_ps) 
)

names(listInput_single_combination_all)[3]<-'combination'
names(listInput_single_combination_all)[4]<-'mofa'

common_mofa_combination=intersect(listInput_single_combination_all$combination, listInput_single_combination_all$mofa)


fname_venn=paste0(out_compare,'ALL_3', 
                  int_params ,'.jpeg')

### def create a venn? 
venn_list=listInput_single_combination_all
create_venn(listInput_single_combination_all, fname_venn, main='MOFA and Single combination and Single ')



params_all=paste0('Single weights: ',
                        paste0(single_weights_Var, collapse=', '  ), '\n n factors: ',length(fns)  
                        )
#### compare MOFA to merged ####
#################################
'cytokine mediated signaling pathway' %in% merged_factors_mofa_sig$Description
'cytokine mediated signaling pathway' %in% single_ps_sig

listInput_single_mofa<-list( mofa=unique(merged_factors_mofa_sig$Description), single=single_ps_sig)
filename = paste0(out_compare,'Venn_MOFA_single_combination', 
                  int_params , cor_t, '.png') 
create_venn(listInput_single_mofa,filename, 
            main=paste0('Mofa and Single combination \n', params_all))
                        

inter_mofa_union<-calculate.overlap(listInput_single_mofa)
unique_single<-unique(inter_mofa_union$a2[!(inter_mofa_union$a2 %in%inter_mofa_union$a3)])
unique_mofa<-unique(inter_mofa_union$a1[!(inter_mofa_union$a1 %in%inter_mofa_union$a3)])

listInput_single_mofa$mofa
###### PLOT UNIQUE ONLY #### 
#####
#####
#####
##
#create_combination_plot(use_mofa=use_mofa, merged_paths_fish=mofa_only_all, combined_p_thresh=0.01, text_add='mofa_uniq')
#create_combination_plot(use_mofa=FALSE, merged_paths_fish=single_only_all, combined_p_thresh=0.01, text_add='single_uniq')

unique(unlist(listInput_all_mods_single), use.names=FALSE)




### eVALUATION ####


retrieve_path_ids <- function(paths_description, go_ids){
  #'
  #' Retrieves the path ids from GO from their descriptions
  #' @param paths_description: path descriptions 
  #' @go_ids key value table with ids and descriptions provided by George Minadakis
  #' @return return the GO ids 
  #paths_description=single_ps_sig
  go_ids_all<-go_ids
  colnames(go_ids_all)<-c('ID', 'Description')
  paths_to_id=data.frame(Description=paths_description)
  # which(go_ids_all$Description=='gene silencing by mirna')
  paths_to_id$Description=standardize_go_names(paths_to_id$Description)
  enrich_rna_single$Description=standardize_go_names(enrich_rna_single$Description)
  enrich_proteins_single$Description=standardize_go_names(enrich_proteins_single$Description)
  
  
  # retrieveGO<-merge(paths_to_id, go_ids_all, by='Description',  all.x=TRUE)
  retrieveGO<-merge(paths_to_id, enrich_rna_single, by='Description',  all.x=TRUE, suffix='rna'); retrieveGO$ID
  retrieveGO<-merge(retrieveGO, enrich_proteins_single, by='Description',  all.x=TRUE, suffix=c('rna','prot')); 
  retrieveGO<-merge(retrieveGO, go_ids_all, by='Description',  all.x=TRUE, suffix=c('','all'));retrieveGO$ID
  
  
  retrieveGO[c('Description','IDrna', 'IDprot', 'ID')] 
  retrieveGO[is.na(retrieveGO$ID),]$ID<-retrieveGO[is.na(retrieveGO$ID),]$IDrna
  retrieveGO[is.na(retrieveGO$ID),]$ID<-retrieveGO[is.na(retrieveGO$ID),]$IDprot
  
  #print(which(is.na(retrieveGO$ID)))
  #retrieveGO[is.na(retrieveGO$ID),]
  
  #paths_to_id$Description[!(paths_to_id$Description %in% enrich_rna_single$Description)]
  
  #enrich_mirnas_single$Description[grep('g1 dna damage checkpoint' ,enrich_mirnas_single$Description)]
  #retrieveGO<-merge(retrieveGO, enrich_proteins, by='Description', all=TRUE)
  
  
  return(retrieveGO$ID)
  
  
  
}



## TRUTHSET #### 
go_ids<-read.csv(paste0(data_dir,'ppmi/ppmi_data/go_pathway_info.txt'), sep='\t')

truth_paths<-read.csv(paste0(data_dir,'ppmi/ppmi_data/known_biomarkers_overlap/GO_Biological_Process_2023_table_gene_malacards.txt' ), sep='\t')
truth_paths_sig<-truth_paths[truth_paths$Adjusted.P.value<0.05,]


### preprocess alll 
truth_paths_sig$Description<-gsub('-', ' ', tolower(truth_paths_sig$Term))
truth_paths_sig$Description= sub( '\\ [^ ]+$', '',truth_paths_sig$Description)
listInput_all_mods_single_l=unique(unlist(listInput_all_mods_single, use.names = FALSE))


listInput_all_mods_single_l_ids<-retrieve_path_ids(listInput_all_mods_single_l, go_ids)
single_ps_sig_ids<-retrieve_path_ids(single_ps_sig , go_ids)


merged_factors_mofa_sig_ids<-retrieve_path_ids(merged_factors_mofa_sig$Description, go_ids)


write.csv(single_ps_sig_ids, paste0(merged_path_file_single, '_GOids.csv'), row.names = FALSE)
write.csv(merged_factors_mofa_sig_ids, paste0(merged_path_file_mofa, '_GOids.csv'), row.names = FALSE)
write.csv(listInput_all_mods_single_l_ids, paste0(path_file_single_union, '_GOids.csv' ), row.names = FALSE)

listInput_single_combination_all$combination
paths_description=listInput_all_mods_single_l

listInput_all_mods_single_l_ids



paths_description=merged_factors_mofa_sig$Description







merged_factors_mofa_sig$Description=gsub('-', ' ', tolower(merged_factors_mofa_sig$Description))
single_ps_sig=gsub('-', ' ', tolower(single_ps_sig))


  
listInput_truth_mofa=append(list(truth_paths_sig$Description),
                            list(merged_factors_mofa_sig$Description)
)
#listInput_truth_mofa=append(listInput_truth_mofa,list(merged_factors_mofa_sig$Description))

# intersect( unlist(listInput_all_mods_single, use.names = FALSE), single_ps_sig)

listInput_truth_mofa=append(listInput_truth_mofa,list(listInput_all_mods_single_l))

listInput_truth_mofa=append(listInput_truth_mofa,list(single_ps_sig))

names(listInput_truth_mofa)=c('malacards', 'mofa', 'union', 'combination' )



fname_venn_all=paste0(out_compare,'VennEuler_Truthset_', 
                      int_params ,'.jpeg')

create_venneuler(listInput_truth_mofa, fname_venn_all)
fname_venn

create_venn(listInput_truth_mofa, fname_venn,
            main=paste('MOFA, union and Malacards \n ',params_all ))


listInput_single_combination

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


unique(listInput_truth_mofa$mofa)

junion=jaccard(unique(listInput_truth_mofa$malacards), unique(listInput_truth_mofa$union))
jmofa=jaccard(unique(listInput_truth_mofa$malacards), unique(listInput_truth_mofa$mofa))
jcomb=jaccard(unique(listInput_truth_mofa$malacards), unique(listInput_truth_mofa$combination))

jaccard(listInput_truth_mofa$malacards, listInput_truth_mofa$combination)

df_j_stats=  format(c( single_weights_Var,adj_weights, junion,jmofa,jcomb ), digits = 2)
df_j_stats
write.table(t(df_j_stats), paste0(out_compare,'jaccard_stats.csv'), append=TRUE, col.names = FALSE)













  ######## THIS IS IF WE ARE USING A DIFFERENT METHOD TO SOURCE PATHS FOR MOFA



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



