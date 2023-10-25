

library(org.Hs.eg.db)
library(edgeR)
source(paste0(script_dir, 'ppmi/utils.R'))
source(paste0(script_dir, 'ppmi/time_utils.R'))

### TODO: run analyze clin vars to load clinvars for later times 




# load this only once..? 
### TODO: ADD PROTEINS TOO!!!! 
process_mirnas=TRUE;
source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se_mirs=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 


process_mirnas=FALSE; source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se_rnas=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 




#### Markers over time:
#### 1. Obtain the markers here either from MOFA OR from deseq 


mode='prognosis'
#mode='prognosis'
## Where to get the molecules from? 
mode_mols='single_time'
model_subtyping<-'MOFA'
mode='diagnosis'
# IN THE DIAGNOSIS MODE we select factors related

#### Markers over time:
#### 1. Obtain the markers here 
fn_sel=4; 
sel_factors_diff
if (mode=='diagnosis'){
  factor=sel_factors[fn_sel]
  sel_factors_mode=sel_factors
}else{
  factor=sel_factors_pd_np3[fn_sel]
  sel_factors_mode=sel_factors_pd_np3
  
}
factor
top_view<-which.max(vars_by_factor[factor,])
top_view

factor=12
factor
top_view<-which.max(vars_by_factor[factor,])
top_view
if (names(top_view)=='miRNA'){
  view='miRNA'; process_mirnas=TRUE; se=se_mirs
  
}else{
  view='RNA'; process_mirnas=FALSE; se=se_rnas 
  
}

#view='miRNA'
se_mirs

#### mofa preprocess
##### Collect molecules that we want to plot over time #### 
#### 1. top MOFA factor molecules
#### 2. top deseq molecules 
#### 3. top timeOmics selected molecules 
mode_mols='MOFA'
if ((mode_mols)=='MOFA'){
  # TODO: function get top x% variables from factor!! 
  f_v<-get_factors(MOFAobject, factors =factor )[[1]]
  ws<-get_weights(MOFAobject, views = view, factors=factor)[[1]]
  cut_high<-0.9; cut_low=1-cut_high
  ws_high<-ws[ws>quantile(ws, cut_high),]
  ws_low<-ws[ws<quantile(ws, cut_low),]
  ws_union<-c(ws_high, ws_low)
  length(ws_union)
  feat_names=names(ws_union)
}else{
  feat_names= sigLRT_genes$gene
  
}






### add clinvars to the requested features too! 
clinvars_to_add<-c('PATNO', 'PATNO_EVENT_ID', 'AGE', 'SEX', 'NHY', 'NP3_TOT', 'COHORT', 'NP3_TOT', 'scopa', 'PDSTATE', 'PD_MED_USE', 
                   'con_putamen')



# create a merged dataframe with all visits to be used downstream. 
# might be better to filer molecules now to save memory..? 
# TODO: which variables are used as id? check when melting feat_names
## only if rna

if (view=='RNA'){
  rownames(se)<-gsub('\\..*', '',rownames(se))
  feat_names<-gsub('\\..*', '',feat_names)
  
}
merged_melt_orig_1<-create_visits_df(se, clinvars_to_add, feat_names = feat_names)
levels(merged_melt_orig_1$variable) # check that the requested variables exist? 





#}

feat_names_ens<-gsub('\\..*', '',feat_names)
feat_names_ens



#feat_names= sigLRT_genes$gene

# Now filter  for the requested molecules 
merged_melt_orig<-merged_melt_orig_1


merged_melt_orig$PATNO
unique(merged_melt_orig$variable)
ens<-gsub('\\..*', '',merged_melt_orig$variable)

#if (view=='RNA'){
#  symb<-get_symbols_vector(ens)
#  merged_melt_orig$symbol<-symb
#  feat_names_ens_ids<-unique(symb)
#}else{
  merged_melt_orig$symbol<-merged_melt_orig$variable
#}


#
### ### NOW match factors to samples
# CREATE GROUPS BY FACTOR 
############################################


sel_cohort=FALSE
# IMPORTANT, IF YOU ADD CONTROLS HERE THEY WILL BE INCLUDED IN THE kmeans grouping!!! 
sel_cohort<-c(1)


if (sel_cohort){
  #'
  #'
  merged_melt=merged_melt_orig[merged_melt_orig$COHORT==sel_cohort, ]
}else{
  merged_melt=merged_melt_orig
}
merged_melt_pd=merged_melt_orig[merged_melt_orig$COHORT==1, ]
merged_melt_ct=merged_melt_orig[merged_melt_orig$COHORT==2, ]




merged_melt_pd<-merged_melt

#### GROUP BY MOFA FACTOR ####
# TODO: FUNCTION
### Wcich mofa run to get factors from??? 
Z <- get_factors(MOFAobject)[[1]]
Z1<-Z[,sel_factors[fn_sel]]
groups_kmeans<-kmeans(Z1, centers=2)
patnos_z1<-gsub('\\_.*', '', names(groups_kmeans$cluster))
groups_kmeans_patnos<-patnos_z1


cluster_by_mofa_factors<-function(MOFAobject, factors,centers=2 ){
  ###
  #' Cluster patients in a mofa object using the specified factors
  #' @param factors which factors to use for the clustering  
  #' @return clusters_by_patno
  #' @
  #' 
  #' 
  
  Z <- get_factors(MOFAobject)[[1]]
  Z1<-Z[,factors]
  groups_kmeans<-kmeans(Z1, centers=centers)
  names(groups_kmeans$cluster)<-gsub('\\_.*', '', names(groups_kmeans$cluster))

  return(groups_kmeans)  
}
factors=factor

factors=factor



groups_from_mofa_factors<-function(merged_melt, MOFAobject, factors ){
  
  #'
  #' @param MOFAobject description
  #'

  ### cluster by one factor 
  
  groups_kmeans<-cluster_by_mofa_factors(MOFAobject, factors=factors, centers=2)
  # OR CHOOSE ALL 

#  Z1_matched<-Z1[match(merged_melt$PATNO,names(groups_kmeans$cluster)) ]
  
  

  high_label<-which.max(groups_kmeans$centers)
  #groups_kmeans$cluster<-ifelse(groups_kmeans$cluster==high_label, 'HighFactor', 'LowFactor')
  pats<-names(groups_kmeans$cluster)
  pats
  ### MATCH the clusters with the dataframe 
  kmeans_matched<-groups_kmeans$cluster[match(merged_melt$PATNO, pats )]
  kmeans_grouping<-factor(kmeans_matched)
  
  groups_from_mofa_factors
  
  
  return(kmeans_grouping)
  
}


### Important: create groups only for the patients.
# cluster patients or controls? 
factors_to_cluster = factor
factors_to_cluster<-sel_factors_diff
factors_to_cluster<-c(sel_factors_np3)

groups_kmeans<-cluster_by_mofa_factors( MOFAobjectPD, factors=factors_to_cluster)
groups_kmeans_diff<-cluster_by_mofa_factors( MOFAobjectPD, factors=sel_factors_diff, centers=2)
groups_kmeans_diff3<-cluster_by_mofa_factors( MOFAobjectPD, factors=sel_factors_diff, centers=3)

merged_melt$kmeans_grouping<-groups_from_mofa_factors(merged_melt, MOFAobjectPD, factors=factors_to_cluster)






na_ps<-unique(merged_melt[!is.na(merged_melt$kmeans_grouping),]$PATNO)
merged_melt_filt<-merged_melt[merged_melt$PATNO %in% na_ps, ]

### EDIT GROUP TO BE PATNO and GROUP!!
# OR JUST GROUP  
merged_melt_filt$VISIT



### Did they change scales???
### Breaks down into two groups based on grouping 



## TODO: plot here ALL the clinical variables by grouping!! 

#### 

merged_melt_filt$VISIT<-as.factor(merged_melt_filt$VISIT)

group_cats<-levels(factor(merged_melt$grouping))








################


### Plot to remove the other group ####
# TAKE THE low group  
# TODO: decide how to take the lowest x and highest x 
### TODO: DO THIS BOTH FOR CONTROLS AND DISEASE ####? 



merged_melt_filt$group<-as.logical(merged_melt_filt$kmeans_grouping)

group_cat='kmeans_grouping'


merged_melt_filt$group<-as.logical(merged_melt_filt[, group_cat])

merged_melt_filt$group<-as.factor(merged_melt_filt[, group_cat] )
group_cats<-levels(merged_melt_filt$group)
## 
# TODO: Function: take a group and 


get_most_sig_over_time<-function(merged_melt_filt_group){
  #''
  #'
  #' @param description
  #'

  merged_melt_filt_group=merged_melt_filt_group[merged_melt_filt_group$VISIT %in% c('BL', 'V08'),]
  merged_melt_filt_group$VISIT<-as.factor(merged_melt_filt_group$VISIT)
  
  
  
  wilcox_stats_group<-merged_melt_filt_group %>%
    group_by(symbol) %>%
    do(w=wilcox.test(value~VISIT, data=.) ) %>%
    summarize(symbol, Wilcox=w$p.value) %>% 
    mutate(p.adj=p.adjust(Wilcox)) %>%
    dplyr::filter(p.adj<0.05) %>%
    arrange(Wilcox, decreasing=FALSE) %>%
    as.data.frame() 
  
  
  most_sig_over_time_group<-wilcox_stats_group[order(wilcox_stats_group$Wilcox),]
  
  
  return(most_sig_over_time_group)
}


merged_melt_ct_two_vis=merged_melt_ct[merged_melt_ct$VISIT %in% c('BL', 'V08'),]
merged_melt_filt_g1=merged_melt_filt[merged_melt_filt$group %in% group_cats[1],]
merged_melt_filt_g2=merged_melt_filt[merged_melt_filt$group %in% group_cats[2],]

most_sig_over_time1<-get_most_sig_over_time(merged_melt_filt_g1)
most_sig_over_time2<-get_most_sig_over_time(merged_melt_filt_g2)
most_sig_over_time_ct<-get_most_sig_over_time(merged_melt_ct)

merged_melt_filt_g2$kmeans_grouping
## TODO: DO IT FOR MULTIPLE GROUPS 
most_sig_over_time<-rbind(most_sig_over_time1, most_sig_over_time2)


######## First find out which of the molecules significantly change over time ####


#### CHOOSE 
merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in%  most_sig_over_time_ct$symbol,]


# they should chgange in pd but not in controls!! 
# remove the ones in the controls
most_sig_over_time<-rbind(most_sig_over_time1, most_sig_over_time2)
most_sig_over_time<-most_sig_over_time%>%
  arrange(Wilcox, decreasing=FALSE)
 
####OUTPUT MOST SIG OVER TIME 


# 
#### 
# TODO: up to here make it into a function 
# 1. Get time variables by group or MULTIPLE GROUPS
# 2. Groups should be by factor OR by multiple factors 
# 3. for now they are by factor!!! 
# 4. save the grouping mode 
factors_to_cluster_s<-paste0(c(factors_to_cluster), collapse='-')

### We get the most significant by group for different factor top variables 
write.csv(most_sig_over_time, paste0(outdir, '/trajectories/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', view,'_',group_cat , '.csv'))

####### CHOOSE 

#most_sig_over_time_deseq = make.names(sigLRT_genes$gene)
#most_sig_over_time_deseq<-make.names(colnames(data.filtered.only.pd))

#most_sig_over_time_deseq



### First answer : CAN THEY DIFFERENTIATE DISEASE CONTROL? 
# TODO: ADD DISEASE CONTROL
merged_melt_ct$kmeans_grouping='CONTROL'


filt_top=TRUE


### PUT THEM ALL TOGETHER IN THE BOXPLOTS 
#merged_melt_all<-rbind(merged_melt_ct, merged_melt_filt_g2_sig)
#merged_melt_all<-rbind(merged_melt_all, merged_melt_filt_g1_sig)




#######################################################
############ TIME TRAJECTORY FOR ALL VISITS ###########
#######################################################



#### Check
mean_data<-merged_melt_filt_g1 %>%
  group_by(grouping, symbol) %>%
  summarise(mean_expr = mean(value, na.rm = TRUE))
mean_data
most_sig_over_time

median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}



merged_melt_ct$kmeans_grouping='CONTROL'
merged_melt_ct$group='CONTROL'
merged_melt_ct$grouping='CONTROL'


merged_melt_filt=rbind(merged_melt_filt,merged_melt_ct )

filt_top=TRUE


if (filt_top){
  
  # TODO: ADD the clinical variables here? 
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% most_sig_over_time$symbol[1:20],]
  
  
  nrow=NULL; height=7
}else{
  merged_melt_filt_most_sig<-merged_melt_filt
  nrow=NULL; height=7
  
}


if (view=='RNA'){
  ens<-as.character(merged_melt_filt_most_sig$symbol)
  symb<-get_symbols_vector(ens)
  merged_melt_filt_most_sig$symbol<-symb
  feat_names_ens_ids<-unique(symb)
}







#### in the boxplots add the groups 
### first controls-- all markers need to be different in controls
### and second in the two groups of disease 
# Boxplots of the grouping too !! 
# TODO: SEPARATE BY PD STATE



ggplot(data = merged_melt_filt_most_sig, aes(x = VISIT, y = value, fill=kmeans_grouping)) + 
  #geom_point(aes(col=VISIT), size = 2) +
  #geom_line(aes(group=PATNO),  col= 'grey') +
  # subgroup should be in the fill parameter!!! 
  geom_boxplot(aes(x=VISIT, fill=kmeans_grouping ))+
  scale_color_viridis_d(option='mako')+
  scale_fill_viridis_d(option='mako')+
  
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +
  
  geom_signif(comparisons = list(c('BL', 'V08')),  
              map_signif_level=TRUE, 
              tip_length = 0, vjust=0.4)+
  
  facet_wrap(. ~ symbol, scales='free_y') +
  
  theme_bw() 
ggsave(paste0(outdir, '/trajectories/boxplots_',factor,'_', view,'_',group_cat,sel_cohort , '.jpeg'), 
       width=12, height=12)


















### BY GROUP ####

# TODO: choose 3 colours grey as control

ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
                                                    fill='group', group='group', colour='group')) + 
  stat_summary(geom = "pointrange", fun.data = median_IQR, 
               position=position_dodge(0))+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1, alpha=0.7) + 
  scale_color_viridis_d(option='magma')+
  facet_wrap(. ~ symbol, scales='free_y', 
             nrow = nrow) +
 # theme_gray()+
  
  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  #
  #geom_signif(comparisons = split(t(combn(levels(merged_melt_filt_most_sig$group), 2)), 
   #                               seq(nrow(t(combn(levels(merged_melt_filt_most_sig$VISIT), 2))))), 
  #            map_signif_level = TRUE, 
  #             tip_length = 0, vjust=0.4)+
  geom_signif(comparisons = list(c('BL', 'V08')), 
              map_signif_level=TRUE, 
              tip_length = 0, vjust=0.3)+
  
  labs(y='logCPM')+
  # legend(legend=c('Low', 'High'))+
  theme(strip.text = element_text(
    size = 10, color = "dark green", face="bold"), 
    axis.title.y =element_text(
      size = 10, color = "dark green", face="bold",), 
    axis.text.x = element_text(
      size = 10 ))



#warnings()
ggsave(paste0(outdir, '/trajectories/trajectory', factor,'_', view,  group_cat,'_',  factors_to_cluster_s, '_', filt_top,sel_cohort,  '.jpeg'), 
       width=7, height=height)




graphics.off()

## color the ones with the highest changes 
# TOP PATIENTS WITH LARGER CHANGES


#top_change<-molecules_change_by_patno[order(molecules_change_by_patno$diff, decreasing = TRUE)[1:20],'PATNO']
#merged_melt_filt_most_sig$TOP=FALSE
#merged_melt_filt_most_sig[merged_melt_filt_most_sig$PATNO %in% top_change,]$TOP<-TRUE
#any(merged_melt_filt_most_sig$TOP)

#### BY PATIENT #####
#p<-ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
#                                                   fill='group', group='group', colour='group')) + 

merged_melt_filt_most_sig$group=merged_melt_filt_most_sig$kmeans_grouping
#p<-ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
#                                                   fill='TOP', group='TOP', colour='TOP')) + 
#
#  geom_point(aes_string(x = 'VISIT', y = 'value', 
#             fill='group', group='group', colour='group' ),size=0.1, alpha=0.5)+
#  geom_line(aes_string(x = 'VISIT', y = 'value', 
#                         group='PATNO', colour='group' ),size=0.1, alpha=0.5)+
#  stat_summary(fun = median, position=position_dodge(width=0), 
#               geom = "line", size = 1) + 
p<-ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
                                                       fill='TOP', group='TOP', colour='TOP')) + 
  
  geom_point(aes_string(x = 'VISIT', y = 'value', 
                        fill='TOP', group='TOP', colour='TOP' ),size=0.1, alpha=0.5)+
  geom_line(aes_string(x = 'VISIT', y = 'value', 
                       group='PATNO', colour='TOP' ),size=0.2, alpha=0.5)+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1) 
### CHANGE OF MOLECULE VS CHANGE OF NP3


add_molecules_changes=FALSE
merged_melt_filt_most_sig$TOP=FALSE

if (add_molecules_changes){
  
  
  levels(merged_melt_filt$symbol)
  
  merged_melt_filt_1  
  ### split by visit 
  molecules_by_visit<-split(merged_melt_filt_1, merged_melt_filt_1$VISIT )
  
  molecules_by_visit2 <- molecules_by_visit %>% 
    imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -PATNO)) %>%
    reduce(full_join, by = "PATNO")
  
  
  molecules_by_visit2
  
  X2=molecules_by_visit2[,paste0('value','_','V08')]
  X1=molecules_by_visit2[,paste0('value','_','BL')]
  
  length(X2)
  length(X1)
  
  
  molecules_by_visit2$log_FC<-(X2-X1)/(X2+X1)
  molecules_by_visit2$diff<-(X2-X1)
  
  #molecules_by_visit2$log_FC<-log2(log(X2)/log(X1))
  
  molecules_by_visit2$log_FC
  
  
  
  molecules_by_visit2
  
  
  
  
  # TOP NEGATIVE CHANGE!
  merged_melt_filt$value
  ### 1. LARGE DIFFERENCES
  # 2. Large changes 
  # 3. large end points 
  top_change<-molecules_change_by_patno[order(molecules_change_by_patno$diff, decreasing = FALSE)[1:20],'PATNO']
  top_change2<-molecules_change_by_patno[order(molecules_change_by_patno$log_FC, decreasing = FALSE)[1:20],'PATNO']
  just_molecules<-merged_melt[merged_melt$variable %in%  c(most_sig_over_time_deseq), ][, c('value', 'PATNO')]
  top_change3<-just_molecules[just_molecules$value< (-0),]$PATNO
  
  
  top_change3
  merged_melt_filt_most_sig$TOP=FALSE
  merged_melt_filt_most_sig[merged_melt_filt_most_sig$PATNO %in% c( top_change2,top_change,top_change3) ,]$TOP<-TRUE
  any(merged_melt_filt_most_sig$TOP)
  
  
  top_molecular_patients<-c( top_change2,top_change,top_change3)
  top_molecular_patients
  
  
}




#### BY PATIENT #####

merged_melt_filt_most_sig$group=merged_melt_filt_most_sig$TOP
merged_melt_filt_most_sig$group=merged_melt_filt_most_sig$kmeans_grouping
merged_melt_filt_most_sig



p<-ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
                                                   fill='group', group='group', colour='group')) + 


  geom_point(aes_string(x = 'VISIT', y = 'value', 
             fill='group', group='group', colour='group' ),size=0.1, alpha=0.5)+
  geom_line(aes_string(x = 'VISIT', y = 'value', 
                         group='PATNO', colour='group' ),size=0.1, alpha=0.5)+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1) + 

  
  theme(legend.position = "none")+
  
  
  
  #scale_color_viridis_d(option='turbo')+
  scale_color_viridis_d(option='magma')+
  
  facet_wrap(. ~ symbol, scales='free_y', 
             nrow = 4) +
  
  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  theme_bw()+ 
  geom_signif(comparisons = list(c('BL', 'V08')), 
              map_signif_level=TRUE, 
              tip_length = 0, vjust=0.4)+
  
  labs(y='logCPM')+
  # legend(legend=c('Low', 'High'))+
  theme(strip.text = element_text(
    size = 13, color = "dark green", face="bold"), 
    axis.title.y =element_text(
      size = 13, color = "dark green", face="bold",), 
    axis.text.x = element_text(
      size = 12 ))

show(p)

#warnings()
ggsave(paste0(outdir, '/trajectories/trajectory_by_pat_', factor,'_', view,  group_cat,'_',  factors_to_cluster_s, '_', filt_top,'_', sel_cohort, '.jpeg'), 
       width=12, height=height)

merged_melt_filt_most_sig 

#df2$change

merged_melt$kmeans_grouping



## collect
# TODO: AUTOMATE this part to collect them all in a function -- retrieve from home? 
factors_to_cluster
sapply(factors_to_cluster )
factor=2
fact2<-as.data.frame(read.csv(paste0(outdir, '/trajectories/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', view,'_',group_cat , '.csv')))

factor=8
fact8<-as.data.frame(read.csv(paste0(outdir, '/trajectories/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', view,'_',group_cat , '.csv')))

factor=9
fact9<-as.data.frame(read.csv(paste0(outdir, '/trajectories/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', view,'_',group_cat , '.csv')))

factor=14
fact14<-as.data.frame(read.csv(paste0(outdir, '/trajectories/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', view,'_',group_cat , '.csv')))



fact2$id <- 2
fact8$id <- 8
fact9$id <- 9
fact14$id <- 14

combined_sig_genes<-rbind(fact2, fact8, fact9, fact14); dim(combined_sig_genes)
combined_sig_genes<-combined_sig_genes[!duplicated(combined_sig_genes$symbol),]

dim(combined_sig_genes)
combined_sig_genes_strict<-combined_sig_genes[combined_sig_genes$p.adj<0.0005,]
#combined_sig_genes_strict<-combined_sig_genes[combined_sig_genes$p.adj<0.0004,]


combined_sig_genes_strict
combined_sig_genes_strict[combined_sig_genes_strict$id==2,]




######## ENRICHMENT BY GENES 

pvalueCutoff<-0.05
gse_2 <- clusterProfiler::enrichGO(fact14$symbol, 
                                   ont=ONT, 
                                   keyType = 'ENSEMBL', 
                                   OrgDb = 'org.Hs.eg.db', 
                                   pvalueCutoff  = pvalueCutoff)



gse_2
results_file_tmp<-paste0(outdir, '/trajectories/enrichment/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', view,'_',group_cat)
enrich_plots<-run_enrichment_plots(gse=gse_2,results_file=results_file_plot, N_DOT=15, N_EMAP=25, text_p=text_p )









### CHANGE OF MOLECULE VS CHANGE OF NP3

### TODO: PLOT FOR THE SAME PATIENTS THE  FUTURE TRAJECTORIES BY GROUPS!! 

#. 1. Add k-means 





sel_visit='V16'
cl_var<-'NP2_TOT'
sel_state = 'OFF'



# TODO FIX 
# this contains all future variable 
df_future_clinvars<-get_future_clinvars(combined_bl_log)


df_to_calc<-get_clinvar_changes(df_future_clinvars, sel_visit = sel_visit,   cl_var=cl_var, sel_state=sel_state)



merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% 'hsa.miR.101.3p',]
merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% 'hsa.miR.101.3p',]
merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% 'hsa.let.7a.3p',]

sel_feature<-most_sig_over_time$symbol[3];sel_feature
#sel_feature<-'ANXA3'
#sel_feature<-'DHRS13'

merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% sel_feature,]





### split by visit 
molecules_by_visit<-split(merged_melt_filt_1, merged_melt_filt_1$VISIT )

molecules_by_visit2 <- molecules_by_visit %>% 
  imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -PATNO)) %>%
  reduce(full_join, by = "PATNO")


molecules_by_visit2

X2=molecules_by_visit2[,paste0('value','_','V08')]
X1=molecules_by_visit2[,paste0('value','_','BL')]

molecules_by_visit2$log_FC<-log(X2/X1)
molecules_by_visit2$diff<-(X2-X1)


scale_change<-df_to_calc[,c( 'diff_scale', 'PATNO',paste0('PDSTATE_', sel_visit ))]
molecules_change_by_patno<-molecules_by_visit2[,c('log_FC','diff', 'PATNO', 'kmeans_grouping_V08')]
molecules_change_by_patno<-merge(molecules_change_by_patno, scale_change, by='PATNO')

hist(scale_change$diff_scale)


### color the top molecular ones too

### WHICH GROUP
kmeans_grouping<-groups_kmeans$cluster
kmeans_grouping<-clusters_patients$cluster
groups_kmeans$centers

names(kmeans_grouping)<-gsub('\\_.*', '', names(kmeans_grouping))

kmeans_grouping=data.frame(kmeans_grouping)
kmeans_grouping$PATNO=rownames(kmeans_grouping)
scale_change$PATNO
kmeans_grouping$PATNO
scale_change_gr<-merge(scale_change, kmeans_grouping, by='PATNO')
scale_change_gr$kmeans_grouping=as.factor(scale_change_gr$kmeans_grouping)
scale_change_gr

ggplot(scale_change_gr, aes(x=diff_scale))+
  geom_histogram(aes(fill=kmeans_grouping))



scale_change<-df_to_calc[,c( 'diff_scale', 'PATNO',paste0('PDSTATE_', sel_visit ))]


molecules_change_by_patno<-molecules_by_visit2[,c('log_FC','diff', 'PATNO', 'kmeans_grouping_V08')]
molecules_change_by_patno<-merge(molecules_change_by_patno, scale_change, by='PATNO')

# Plot the absolute difference between
# Diff

molecules_change_by_patno[, paste0('PDSTATE_', sel_visit )]<-factor(molecules_change_by_patno[, paste0('PDSTATE_', sel_visit )])
colnames(molecules_change_by_patno)
ggplot(molecules_change_by_patno[molecules_change_by_patno$kmeans_grouping_V08!='CONTROL',], 
       aes(x=log_FC, y=diff_scale))+
  geom_point(aes(color=kmeans_grouping_V08))+
  geom_smooth(method = "lm")+
  facet_wrap(as.formula(paste0('~ PDSTATE_', sel_visit)), nrow=3)+
  labs(title=paste(sel_feature, cl_var))


ggsave(paste0(outdir, '/trajectories/change/change_', factor, '_',sel_feature,'_', cl_var, '_',sel_visit,sel_state,'.jpeg'), 
       width=6, height=5)

molecules_change_by_patno























