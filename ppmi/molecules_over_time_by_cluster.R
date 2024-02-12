

library(org.Hs.eg.db)
library(edgeR)
source(paste0(script_dir, 'ppmi/utils.R'))
source(paste0(script_dir, 'ppmi/time_utils.R')); 
source(paste0(script_dir, 'ppmi/plotting_utils.R'))
source(paste0(script_dir, 'ppmi/mofa_utils.R'))

#source(paste0(script_dir, 'ppmi/cluster_comparisons.R'))

#get_de_features_by_group
#library('EnsDb.Hsapiens.v79')
### TODO: run analyze clin vars to load clinvars for later times 


factors_to_cluster_s<-paste0(c(which(all_fs_diff[,y_clust])), collapse='-')
factors_to_cluster_s

# load this only once..? 
process_mirnas=TRUE;
source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se_mirs=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 


process_mirnas=FALSE; source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
#se_rnas=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 
#se_rnas=se_filt_corrected; # load the data after batch correction 
# TODO: Correct also for neutrophils to see the trajectories?
vsd_cor_l=loadRDS(vst_cor_all_vis)
vsd_cor_filt<-filter_se(vsd_cor_l, VISIT = c('BL','V04', 'V06', 'V08'), sel_coh = sel_coh, sel_sub_coh = sel_subcoh)

se_rnas<-vsd_cor_filt

# TODO: load proteins 
VISIT=c('BL', 'V04', 'V06', 'V08')
process_mirnas=FALSE; source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
prot_vsn_se_filt_file
datalist<-loadRDS(prot_vsn_se_filt_file)
vsn_mat_proteins<-datalist[[1]]
se_filt_proteins<-datalist[[2]]





#### Markers over time:
#### 1. Obtain the markers here either from MOFA OR from deseq 


mode='prognosis'
#mode='prognosis'
## Where to get the molecules from? 
mode_mols='single_time'
model_subtyping<-'MOFA'
#mode='diagnosis'
# IN THE DIAGNOSIS MODE we select factors related

#### Markers over time:
#### 1. Obtain the markers here 

fn_sel=23; 
if (mode=='diagnosis'){
  factor=sel_factors[fn_sel]
  sel_factors_mode=sel_factors
}else{
  factor=sel_factors_pd_np3[fn_sel]
  sel_factors_mode=sel_factors_pd_np3
  
}


factor=23
top_view<-which.max(vars_by_factor[factor,])


names(top_view)<-'miRNA'; keep_all_feats=TRUE
names(top_view)<-'RNA'; keep_all_feats=TRUE 
names(top_view)<-'proteomics_plasma'; keep_all_feats=TRUE 

top_factor_feats_rna<-select_top_bottom_perc(MOFAobject=MOFAobject,view='RNA', factor, top_fr = 0.01)
top_factor_feats_mirna<-select_top_bottom_perc(MOFAobject=MOFAobject,view='miRNA', factor, top_fr = 0.01)
top_factor_feats_pl<-select_top_bottom_perc(MOFAobject=MOFAobject,view='proteomics_plasma', factor, top_fr = 0.01)

top_factor_feats_pl

if (names(top_view)=='miRNA'){
  view='miRNA'; process_mirnas=TRUE; se=se_mirs
 top_factor_feats= top_factor_feats_mirna
}else if (names(top_view)=='RNA') {
  view='RNA'; process_mirnas=FALSE; se=se_rnas 
    top_factor_feats = top_factor_feats_rna

}else if (names(top_view)=='proteomics_plasma') {
  view='proteomics_plasma';se=se_filt_proteins
  top_factor_feats = top_factor_feats_pl
  
  }

se_filt_proteins

#### mofa preprocess
##### Collect molecules that we want to plot over time #### 
#### 1. top MOFA factor molecules
#### 2. top deseq molecules 
#### 3. top timeOmics selected molecules 

  # TODO: function get top x% variables from factor!! 
  





feats_rna_all<-rownames(get_weights(MOFAobject, views = 'RNA', factors=factor)[[1]])
feats_mirna_all<-rownames(get_weights(MOFAobject, views = 'miRNA', factors=factor)[[1]])
feats_pl_all<-rownames(get_weights(MOFAobject, views = 'proteomics_plasma', factors=factor)[[1]])

rownames(se_rnas)<-gsub('\\..*', '',rownames(se_rnas))
feats_rna_all<-gsub('\\..*', '',feats_rna_all)
top_factor_feats_rna<-gsub('\\..*', '',top_factor_feats_rna)

### add clinvars to the requested features too! 
clinvars_to_add<-c('PATNO', 'PATNO_EVENT_ID', 'AGE', 'SEX', 'NHY','NP2PTOT', 'NP3TOT', 'COHORT', 'scopa', 'PDSTATE', 'PD_MED_USE', 
                   'con_putamen')



# create a merged dataframe with all visits to be used downstream. 
# might be better to filer molecules now to save memory..? 
# TODO: which variables are used as id? check when melting feat_names
## only if rna

if (view=='proteomics_plasma' || view=='proteomics_csf' ){
  feat_names<-gsub('_proteomics.*', '',feat_names)
  top_factor_feats<-gsub('_proteomics.*', '',top_factor_feats)
  feats_pl_all=gsub('_proteomics.*', '',feats_pl_all)
}



#rescale_option
merged_melt_filt_rna<-create_merged_melt(se_rnas, feats_rna_all, view='RNA', MOFAobject_clusts=MOFAobject_sel)
merged_melt_filt_pl<-create_merged_melt(se_filt_proteins, feats_pl_all, view='proteomics_plasma', MOFAobject_clusts=MOFAobject_sel, 
                run_cpm = FALSE)

#merged_melt_filt_mirna<-create_merged_melt(se_mirnas, feats_mirna_all, view='miRNA')


### List of merged per group
merged_melt_filt=merged_melt_filt_pl
merged_melt_filt$symbol


merged_melt_groups<-merged_melt_filt %>% 
  split(~group)


merged_melt_groups

library(stringr)
# load DE genes from cluster_comparisons file
# load also time related features?? 
# TODO: load these from the other file

deseq_all<- vector("list", length = 3)

de_files<-paste0(deseq_params, '/', prefix, 'de_cluster_' )

#de_file<-paste0(cluster_params_dir, '/',prefix, 'de_cluster_', cluster_id , '.csv')

for (cluster_id in 1:3){
    deseq2ResDF<-read.csv(paste0(de_files,  cluster_id , '.csv'), row.names=1 )
    deseq_all[[cluster_id]]<-deseq2ResDF[deseq2ResDF$mofa_sign %in% 'Significant',]
   # print(head(deseq_all[[cluster_id]]))
    deseq_all_names <- lapply(deseq_all, function(x){
      return(  gsub('\\..*', '',rownames(x))   )  })
    names(deseq_all_names)<-paste0('SG', 1:length(deseq_all_names))
}
names(deseq_all_names)
# check that they are in mofa

sel_feats<-unique(merged_melt_filt$symbol)
sel_feats
intersect(deseq_all_names$SG2, sel_feats)
deseq_all_names$SG1 %in% sel_feats
deseq_all_names
de_group_vs_control<-lapply(deseq_all_names, function(deseq_clust){ deseq_clust[deseq_clust %in% sel_feats ] }) # filter de group vs control for selected
de_group_vs_control
## significant and varying with time 
de_group_vs_control[[1]]


most_sig_over_time_list<-lapply(1:3, function(clust_id){
  if (length(de_group_vs_control[[clust_id]])>0){
   get_most_sig_over_time( merged_melt_groups[[clust_id]][merged_melt_groups[[clust_id]]$symbol %in%
    de_group_vs_control[[clust_id]] ,] )
}}

)

# they should chgange in pd but not in controls!! 
# remove the ones in the controls
most_sig_over_time<-do.call(rbind, most_sig_over_time_list)


## ECHECK DIRECTORIIONS
#geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= c('CSAD'), keytype = c("SYMBOL"), columns=c('SYMBOL', 'GENEID'))






## CHOOSE MOLECULES HERE ####
####OUTPUT MOST SIG OVER TIME 
#most_sig_over_time<-most_sig_over_time[!(most_sig_over_time$symbol %in% most_sig_over_time_ct$symbol),]

#de_group_vs_control_and_time<-intersect(most_sig_over_time$symbol, de_group_vs_control$symbol)

#### Choose molecules DE in group 1, but most sig over time in all? 


choose_group<-1

#most_sig_over_time1$symbol

de_group_vs_control_and_time<-intersect(most_sig_over_time_list[[choose_group]]$symbol, de_group_vs_control[[choose_group]])



#'ENSG00000205683' %in% de_group_vs_control_and_time
##### CHECK MONOTONICITY AND REMOVE ####

#### WHATEVER SET WE CHOOSE WE REMOVE the non monotonic genes 

# 
#### 
# TODO: up to here make it into a function 
# 1. Get time variables by group or MULTIPLE GROUPS
# 2. Groups should be by factor OR by multiple factors 
# 3. for now they are by factor!!! 
# 4. save the grouping mode 
factors_to_cluster_s<-paste0(c(which(all_fs_diff[,y_clust])), collapse='-')
factors_to_cluster_s



if (view=='RNA'){
  ens<-as.character(most_sig_over_time$symbol)
  symb<-get_symbols_vector(ens)
  most_sig_over_time$GENE_SYMBOL<-symb
  #feat_names_ens_ids<-unique(symb)
  
  
  
}


choose_group=2
### We get the most significant by group for different factor top variables 
write.csv(most_sig_over_time1, paste0(outdir, '/trajectories/most_sig_over_time_',factor,'feats_', keep_all_feats, '_cl_fs_',factors_to_cluster_s,'_', 
                                     view,'_' ,
                                    'de_group_',  choose_group, 
                                     '.csv'))

dim(de_group_vs_control1)
write.csv(de_group_vs_control1, paste0(outdir, '/trajectories/most_sig_vs_control_',factor,'feats_', keep_all_feats, '_cl_fs_',factors_to_cluster_s,'_', 
                                      view,'_' ,
                                      'de_group_',  choose_group, 
                                      '.csv'))



### PUT THEM ALL TOGETHER IN THE BOXPLOTS 



#######################################################
############ TIME TRAJECTORY FOR ALL VISITS ###########
#######################################################
merged_melt_filt = merged_melt_filt_pl

filt_top='top'; filt_top='selected'; filt_top='all' 
filt_top='top'; 

#selected_rnas<-unique(mirna_targets_edgel$Subcategory )# top10 factor 2)

selected_rnas<-top_factor_feats_rna
if (view=='miRNA'){
  selected_feats<-selected_mirs
}else if ((view=='RNA')){
  selected_feats<-selected_rnas
}

if (view=='RNA'){
  ens<-as.character(merged_melt_filt$symbol)
  symb<-get_symbols_vector(ens)
  merged_melt_filt$GENE_SYMBOL<-symb
}else{
  merged_melt_filt$GENE_SYMBOL<-merged_melt_filt$symbol
  
}


# unique in 1? 

de_group_vs_control1_unique_top<-intersect(top_factor_feats_rna,deseq_all_top$SG1[!(deseq_all_top$SG1 %in% c(deseq_all_top$SG2, deseq_all_top$SG3) )] )
de_group_vs_control2_unique_top<-intersect(top_factor_feats_rna,deseq_all_top$SG2[!(deseq_all_top$SG2 %in% c(deseq_all_top$SG1, deseq_all_top$SG3) )] )
de_group_vs_control3_unique_top<-intersect(top_factor_feats_rna,deseq_all_top$SG3[!(deseq_all_top$SG3 %in% c(deseq_all_top$SG1, deseq_all_top$SG2) )] )




### Rules: 
# 1. DE and time difference 
# 2.  De only
# 3. de in specific group unique
# 4. de and in selected pathways 
head(deseq_all_names$SG2)


choose_group=1

merged_melt_filt$symbol
merged_melt_filt_most_sig$symbol
#top_factor_feats = top_factor_feats_pl
merged_melt_filt
merged_melt_filt$symbol
top_factor_feats
top_factor_feats



if (filt_top=='top'){

  # TODO: ADD the clinical variables here? 
#  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$GENE_SYMBOL %in% top10_selected_paths,]
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% de_group_vs_control2[1:30],]
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% top_factor_feats,]
  
  nrow=NULL; height=7; width=18
  nrow=6
}else if (filt_top=='selected'){
  ## keep specific mirs for presentation 
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% selected_feats,]
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$GENE_SYMBOL %in% selected_feats,]
  
  nrow=1; height=3; width=9;
  
}else{
  merged_melt_filt_most_sig<-merged_melt_filt
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% de_group_vs_control_and_time2,]
  #merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% de_group_vs_control_and_time2,]
  
  nrow=NULL; height=12; 
}
merged_melt_filt_most_sig$symbol
unique(merged_melt_filt_most_sig$GENE_SYMBOL)
view
if (view=='RNA'){
  ens<-as.character(merged_melt_filt_most_sig$symbol)
  symb<-get_symbols_vector(ens)
  merged_melt_filt_most_sig$symbol<-symb
  #feat_names_ens_ids<-unique(symb)
  
  
  
}


#### in the boxplots add the groups 
### first controls-- all markers need to be different in controls
### and second in the two groups of disease 
# Boxplots of the grouping too !! 
# TODO: SEPARATE BY PD STATE

#plot_molecular_trajectories(merged_melt_filt_most_sig)



median_IQR <- function(x) {
  data.frame(y = median(x, na.rm = TRUE), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}

merged_melt_filt_most_sig$month <-as.factor(as.numeric( mapvalues(as.character(merged_melt_filt_most_sig$VISIT), 
                                                   from= names(EVENT_MAP), 
                                                   to=unlist(EVENT_MAP, use.names=FALSE))))

factors_to_clust_s<-paste0(unlist(factors_to_clust), collapse='_')
trajectory_fname<-paste0(outdir, '/trajectories/trajectory', factor,'_',keep_all_feats,'_', view, 
       '_',  factors_to_clust_s, filt_top,
       'cluster_',choose_group,'.jpeg')

trajectory_fname
plot_molecular_trajectories_line(merged_melt_filt_most_sig,x='month', trajectory_fname = trajectory_fname )



length(feats_rna_all)
#top_factor_feats_rna

























