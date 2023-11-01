

library(org.Hs.eg.db)
library(edgeR)
source(paste0(script_dir, 'ppmi/utils.R'))
source(paste0(script_dir, 'ppmi/time_utils.R'))
source(paste0(script_dir, 'ppmi/plotting_utils.R'))

library('EnsDb.Hsapiens.v79')
### TODO: run analyze clin vars to load clinvars for later times 




# load this only once..? 
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
if (mode=='diagnosis'){
  factor=sel_factors[fn_sel]
  sel_factors_mode=sel_factors
}else{
  factor=sel_factors_pd_np3[fn_sel]
  sel_factors_mode=sel_factors_pd_np3
  
}


factor=12
top_view<-which.max(vars_by_factor[factor,])
#top_view='miRNA'


#names(top_view)<-'RNA'
if (names(top_view)=='miRNA'){
  view='miRNA'; process_mirnas=TRUE; se=se_mirs
  
}else{
  view='RNA'; process_mirnas=FALSE; se=se_rnas 
  
}

view
se_mirs

#### mofa preprocess
##### Collect molecules that we want to plot over time #### 
#### 1. top MOFA factor molecules
#### 2. top deseq molecules 
#### 3. top timeOmics selected molecules 
mode_mols='MOFA'
keep_all_feats=TRUE
keep_all_feats=FALSE

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
  if (keep_all_feats){
    feat_names=rownames(ws)
    
  }
  
}else{
  feat_names= sigLRT_genes$gene
  
}

### if use all 
#feat_names=rownames(se)

length(feat_names)

### add clinvars to the requested features too! 
clinvars_to_add<-c('PATNO', 'PATNO_EVENT_ID', 'AGE', 'SEX', 'NHY','NP2PTOT', 'NP3TOT', 'COHORT', 'scopa', 'PDSTATE', 'PD_MED_USE', 
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
merged_melt_orig$symbol<-merged_melt_orig$variable



#
### ### NOW match factors to samples
# CREATE GROUPS BY FACTOR 
############################################



# IMPORTANT, IF YOU ADD CONTROLS HERE THEY WILL BE INCLUDED IN THE kmeans grouping!!! 
sel_cohort<-c(1)
  sel_cohort=FALSE

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
y_clust='NP2PTOT'
clust_name= paste0(y_clust, '_clust')

names(all_clusts_mofa)
MOFAobject@samples_metadata
MOFAobject_clusts=MOFAobjectPD
groups_from_mofa_factors<-function(patnos, MOFAobject_clusts, y_clust ){
  
  #' Obtain the molecular clusters from the mofa object 
  #' 
  #' @param MOFAobject description
  #' @param

  ### cluster by one factor 
  clust_name= paste0(y_clust, '_clust')
  sm_pd<-MOFAobject_clusts@samples_metadata;
  groups_kmeans<-sm_pd[, clust_name]
  pats<-sm_pd$PATNO
  kmeans_matched<-groups_kmeans[match(patnos, pats )]
  kmeans_grouping<-factor(kmeans_matched)
  kmeans_grouping
  return(kmeans_grouping)
  
}


merged_melt$kmeans_grouping<-groups_from_mofa_factors(patnos=merged_melt$PATNO, MOFAobject_clusts= MOFAobjectPD, y_clust)
merged_melt$kmeans_grouping=as.numeric(merged_melt$kmeans_grouping)
merged_melt[merged_melt$COHORT%in%c(2), 'kmeans_grouping']<-'HC'

# ADD LABELS FOR controls
merged_melt$kmeans_grouping<-as.factor(merged_melt$kmeans_grouping)

na_ps<-unique(merged_melt[!is.na(merged_melt$kmeans_grouping),]$PATNO)
merged_melt_filt<-merged_melt[merged_melt$PATNO %in% na_ps, ]
merged_melt_filt$VISIT<-as.factor(merged_melt_filt$VISIT)



################


### Plot to remove the other group ####
# TAKE THE low group  
# TODO: decide how to take the lowest x and highest x 
### TODO: DO THIS BOTH FOR CONTROLS AND DISEASE ####? 



merged_melt_filt$group<-merged_melt_filt$kmeans_grouping

group_cat='kmeans_grouping'
merged_melt_filt$group<-as.factor(merged_melt_filt[, group_cat] )
group_cats<-levels(merged_melt_filt$group)
## 
# TODO: Function: take a group and DE it with controls 




### List of merged per group

merged_melt_groups<-merged_melt_filt %>% 
  split(~group)
  


merged_melt_filt_g1=merged_melt_groups[[1]]
merged_melt_filt_g1_ct<-rbind(merged_melt_groups[[1]], merged_melt_groups[['HC']])
merged_melt_filt_g2=merged_melt_groups[[2]]
merged_melt_filt_g2_ct=rbind(merged_melt_groups[[2]], merged_melt_groups[['HC']])
merged_melt_filt_g3=merged_melt_groups[[3]]
merged_melt_filt_g3_ct=rbind(merged_melt_groups[[3]], merged_melt_groups[['HC']])



library(stringr)


de_group_vs_control1<-get_de_features_by_group(merged_melt_filt_g1_ct, var_to_diff='kmeans_grouping') %>%  dplyr::filter(str_detect(symbol, "^ENS|^hsa"))
de_group_vs_control2=NULL
de_group_vs_control2<-get_de_features_by_group(merged_melt_filt_g2_ct, var_to_diff='kmeans_grouping')%>% dplyr::filter(str_detect(symbol, "^ENS|^hsa"))
de_group_vs_control3<-get_de_features_by_group(merged_melt_filt_g3_ct, var_to_diff='kmeans_grouping')%>% dplyr::filter(str_detect(symbol, "^ENS|^hsa"))


### CHECK MONOTONICITY FOR common de genes
# idea: if there is a change: 
# then check that the change is  in the same direction 

de_merged<-merge(de_group_vs_control1,de_group_vs_control3, by='symbol')
de_merged2<-merge(de_group_vs_control1,de_group_vs_control2, by='symbol' )
de_merged_to_remove1<-de_merged[!( de_merged$direction.x ==de_merged$direction.y), ]; 
de_merged_to_remove2<-de_merged2[!(de_merged2$direction.x ==de_merged2$direction.y), ]

dim(de_merged)
dim(de_merged_to_remove)
dim(de_group_vs_control1)

de_merged_to_remove=c(de_merged_to_remove1$symbol,de_merged_to_remove2$symbol )
de_merged_to_remove




inter_de<-intersect(de_group_vs_control1$symbol,de_group_vs_control2$symbol )
union_de<-list(group1=de_group_vs_control1$symbol, group2=de_group_vs_control2$symbol, group3=de_group_vs_control3$symbol)
de_group_vs_control<-rbind(de_group_vs_control1, de_group_vs_control2, de_group_vs_control3 )
de_group_vs_control_inter<-intersect(de_group_vs_control1$symbol, de_group_vs_control2$symbol)

#### Venn diagram 
nclusts<-max(unique(group_cats[group_cats!='HC']))
fname_venn=paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/','venn_de_per_group_',keep_all_feats,'.png')
create_venn(venn_list = union_de, fname_venn =fname_venn,main =paste0( ' DE molecules for each molecular cluster' ))


#### GET DE over time from the de features 
length(unique(de_group_vs_control1$symbol))
length(unique(de_group_vs_control2$symbol))
length(unique(de_group_vs_control3$symbol))


most_sig_over_time1<-get_most_sig_over_time(merged_melt_filt_g1[merged_melt_filt_g1$symbol %in% de_group_vs_control1$symbol ,] %>% 
                                            dplyr::filter(!(symbol %in% de_merged_to_remove) ))
most_sig_over_time2=NULL
most_sig_over_time2$symbol=NULL
most_sig_over_time2<-get_most_sig_over_time(merged_melt_filt_g2[merged_melt_filt_g2$symbol %in% de_group_vs_control2$symbol,] %>% 
                                              dplyr::filter(!(symbol %in% de_merged_to_remove) ) )
most_sig_over_time3<-get_most_sig_over_time(merged_melt_filt_g3[merged_melt_filt_g3$symbol %in% de_group_vs_control3$symbol,]%>% 
                                              dplyr::filter(!(symbol %in% de_merged_to_remove) )
                                            )

most_sig_over_time1

length(most_sig_over_time1$symbol); length(most_sig_over_time2$symbol); length(most_sig_over_time3$symbol)

union_de_time<-list(group1=most_sig_over_time1$symbol, group2=most_sig_over_time2$symbol, group3=most_sig_over_time3$symbol)
fname_venn=paste0(outdir, '/clustering/', clust_name, '/',nclusts,'/', rescale_option, '/','venn_de_per_group_time', keep_all_feats,'.png')
create_venn(venn_list = union_de_time, fname_venn =fname_venn,main =paste0( ' DE molecules with temporal and monotonic change
                                                                            between the clusters and controls' ))




intersect(most_sig_over_time1$symbol, de_group_vs_control1$symbol)

most_sig_over_time_ct<-get_most_sig_over_time(merged_melt_ct)

## TODO: DO IT FOR MULTIPLE GROUPS 
most_sig_over_time<-rbind(most_sig_over_time1, most_sig_over_time2)



######## First find out which of the molecules significantly change over time ####


#### CHOOSE 
merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in%  most_sig_over_time2$symbol,]


# they should chgange in pd but not in controls!! 
# remove the ones in the controls
most_sig_over_time<-rbind(most_sig_over_time1, most_sig_over_time2)
most_sig_over_time<-most_sig_over_time%>%
  arrange(Wilcox, decreasing=FALSE)

## ECHECK DIRECTORIIONS
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= c('CSAD'), keytype = c("SYMBOL"), columns=c('SYMBOL', 'GENEID'))

de_group_vs_control3[de_group_vs_control3$symbol==geneIDs1$GENEID,]
de_group_vs_control2[de_group_vs_control2$symbol==geneIDs1$GENEID,]

de_group_vs_control1[de_group_vs_control1$symbol==geneIDs1$GENEID,]





## CHOOSE MOLECULES HERE ####
####OUTPUT MOST SIG OVER TIME 
most_sig_over_time<-most_sig_over_time[!(most_sig_over_time$symbol %in% most_sig_over_time_ct$symbol),]

de_group_vs_control_and_time<-intersect(most_sig_over_time$symbol, de_group_vs_control$symbol)

#### Choose molecules DE in group 1, but most sig over time in all? 

most_sig_over_time1

choose_group<-1
if (choose_group==1){
  de_group_vs_control_and_time<-intersect(most_sig_over_time1$symbol, de_group_vs_control1$symbol)
  
  
}else if(choose_group==2){
  de_group_vs_control_and_time<-intersect(most_sig_over_time2$symbol, de_group_vs_control2$symbol)
  
}else if(choose_group==3){
  de_group_vs_control_and_time<-intersect(most_sig_over_time3$symbol, de_group_vs_control3$symbol)
  
}

'ENSG00000205683' %in% de_group_vs_control_and_time
##### CHECK MONOTONICITY AND REMOVE ####

#### WHATEVER SET WE CHOOSE WE REMOVE the non monotonic genes 
de_group_vs_control_and_time2<-de_group_vs_control_and_time[!(de_group_vs_control_and_time%in% de_merged_to_remove)]


length(de_group_vs_control_and_time2)


# 
#### 
# TODO: up to here make it into a function 
# 1. Get time variables by group or MULTIPLE GROUPS
# 2. Groups should be by factor OR by multiple factors 
# 3. for now they are by factor!!! 
# 4. save the grouping mode 
factors_to_cluster_s<-paste0(c(which(all_fs_diff[,y_clust])), collapse='-')
factors_to_cluster_s
### We get the most significant by group for different factor top variables 
write.csv(most_sig_over_time, paste0(outdir, '/trajectories/most_sig_over_time_',factor, '_cl_fs_',factors_to_cluster_s,'_', 
                                     view,'_',group_cat ,
                                    'de_group_',  choose_group, 
                                     '.csv'))



### PUT THEM ALL TOGETHER IN THE BOXPLOTS 
#merged_melt_all<-rbind(merged_melt_ct, merged_melt_filt_g2_sig)
#merged_melt_all<-rbind(merged_melt_all, merged_melt_filt_g1_sig)




#######################################################
############ TIME TRAJECTORY FOR ALL VISITS ###########
#######################################################



de_group_vs_control_and_time2

filt_top=TRUE
de_group_vs_control_and_time

unique(de_group_vs_control_and_time2)

if (filt_top){

  # TODO: ADD the clinical variables here? 
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% most_sig_over_time$symbol[1:20],]
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% de_group_vs_control_and_time2[1:25],]
  
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% de_group_vs_control_and_time2,]
  
  nrow=NULL; height=7
}else{
  merged_melt_filt_most_sig<-merged_melt_filt
  nrow=NULL; height=7
  
}

de_group_vs_control_and_time2
if (view=='RNA'){
  ens<-as.character(merged_melt_filt_most_sig$symbol)
  symb<-get_symbols_vector(ens)
  merged_melt_filt_most_sig$symbol<-symb
  feat_names_ens_ids<-unique(symb)
}


#merged_melt_filt_most_sig_g1<-

print(merged_melt_filt_most_sig$NP2PTOT )


#### in the boxplots add the groups 
### first controls-- all markers need to be different in controls
### and second in the two groups of disease 
# Boxplots of the grouping too !! 
# TODO: SEPARATE BY PD STATE

plot_molecular_trajectories(merged_melt_filt_most_sig)



median_IQR <- function(x) {
  data.frame(y = median(x), # Median
             ymin = quantile(x)[2], # 1st quartile
             ymax = quantile(x)[4])  # 3rd quartile
}




# TODO: choose 3 colours grey as control
add_patient_lines=FALSE
merged_melt_filt_most_sig<-merged_melt_filt_most_sig %>% dplyr::filter(VISIT%in% c('BL', 'V04', 'V08'))
p<-ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
                                                       fill='group', group='group', colour='group')) 
if (add_patient_lines){
  p<- p+geom_line(aes_string(x = 'VISIT', y = 'value', 
                             group='PATNO', colour='group' ),size=0.1, alpha=0.5)
}


p=p+ stat_summary(geom = "pointrange", fun.data = median_IQR, 
                  position=position_dodge(0))+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1, alpha=0.7) + 
  scale_color_viridis_d(option='magma')+
  facet_wrap(. ~ symbol, scales='free_y', 
             nrow = nrow) +
  
  #geom_signif(comparisons = list(c('BL', 'V08')), 
  #            map_signif_level=TRUE, 
  #            tip_length = 0, vjust=0.3)+
  
  labs(y='logCPM')+
  # legend(legend=c('Low', 'High'))+
  theme(strip.text = element_text(
    size = 10, color = "dark green", face="bold"), 
    axis.title.y =element_text(
      size = 10, color = "dark green", face="bold",), 
    axis.text.x = element_text(
      size = 10 ))



#warnings()
ggsave(paste0(outdir, '/trajectories/trajectory', factor,'_',keep_all_feats,'_', view,  group_cat,'_',  factors_to_cluster_s, '_', filt_top,sel_cohort,
              'cluster_',choose_group,'.jpeg'), 
       width=20, height=height)




graphics.off()

















merged_melt_filt_most_sig$value=as.numeric(merged_melt_filt_most_sig$value)
merged_melt_filt_most_sig$NP2PTOT=as.numeric(merged_melt_filt_most_sig$NP2PTOT)

merged_melt_filt_most_sig_g1<-merged_melt_filt_most_sig[merged_melt_filt_most_sig$kmeans_grouping==1,]


merged_melt_filt_most_sig_g1$NP2PTOT
crtest<-cor.test(merged_melt_filt_most_sig_g1$value, merged_melt_filt_most_sig_g1$NP2PTOT)
var1<-merged_melt_filt_most_sig_g1 %>%
  dplyr::filter(variable == 'ENSG00000149527')

all_cors_with_scale<-merged_melt_filt_most_sig_g1%>%
  group_by(variable) %>%
  mutate(value=as.numeric(value))%>%
  mutate(NP2PTOT=as.numeric(NP2PTOT))%>%
  dplyr::summarize(cor= cor.test(value, NP2PTOT, use = "complete.obs")$estimate)%>%
  arrange(desc(cor))

all_cors_with_scale

var1$NP2PTOT_cut<-cut(var1$NP2PTOT,breaks=4)
ggplot(var1,aes( x=NP2PTOT_cut, y=value))+
  geom_boxplot()











### collect all results 
sig_genes_all_factors=data.frame()
factors_to_fetch<-c(sel_factors, sel_factors_pd_np3)
factors_to_fetch<-sel_factors_diff

view='RNA'
factors_to_fetch=c(1,4,6,8,9,11,12)
all_files<-sapply(factors_to_fetch, function(factor){
  ### collect all moelcules for each factor 
  top_view<-which.max(vars_by_factor[factor,])
  view=names(top_view)
  
  most_sig_file1<-paste0(outdir, '/trajectories/most_sig_over_time_',factor,'_', view,'_',group_cat , '.csv')
  return(most_sig_file1)
  
})
all_files[[3]]

file.exists(all_files[[3]])

all_sig_genes<-lapply(all_files, function(file){
  if (file.exists(file)){
    print('found_file')
    sig_genes_f<-as.data.frame(read.csv(file))
    return(sig_genes_f)
  }
})

  

all_sig_genes2 <- do.call("rbind", all_sig_genes)
all_sig_genes2<-all_sig_genes2%>%
  arrange(Wilcox, decreasing=FALSE)


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










