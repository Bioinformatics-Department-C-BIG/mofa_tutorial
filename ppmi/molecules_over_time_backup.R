
library(org.Hs.eg.db)
library(edgeR)
source(paste0(script_dir, 'ppmi/utils.R'))

### TODO: run analyze clin vars to load clinvars for later times 



#### Markers over time:
#### 1. Obtain the markers here either from MOFA OR from deseq 


mode='diagnosis'
## Where to get the molecules from? 
mode_mols='single_time'
model_subtyping<-'MOFA'
# IN THE DIAGNOSIS MODE we select factors related

#### Markers over time:
#### 1. Obtain the markers here 
fn_sel=1; 
if (mode=='diagnosis'){
  factor=sel_factors[fn_sel]
  
}else{
  factor=sel_factors_pd_np3[fn_sel]
}


factor
top_view<-which.max(vars_by_factor[factor,])
top_view
if ((top_view)=='miRNA'){
  view='miRNA'; process_mirnas=TRUE
  
}else{
  view='RNA'; process_mirnas=FALSE
  
}


source(paste0(script_dir, 'ppmi/config.R'));deseq_file;
se=load_se_all_visits(input_file = input_file, combined=combined_bl_log); 

create_visits_df<-function(se){
  se_filt_V08<-filter_se(se, VISIT='V08', sel_coh,sel_ps)
  se_filt_BL<-filter_se(se, VISIT='BL', sel_coh,sel_ps)
  se_filt_V06<-filter_se(se, VISIT='V06', sel_coh,sel_ps)
  se_filt_V04<-filter_se(se, VISIT='V04', sel_coh,sel_ps)
  z
}


#if (view==proteomics){
#  
#}
### IF proteins??? 


levels(merged_melt_orig_1$variable)


#### mofa preprocess
##### Collect molecules that we want to plot over time #### 
#### 1. top MOFA factor molecules
#### 2. top deseq molecules 
#### 3. top timeOmics selected molecules 
mode_mols='MOFA'
if ((mode_mols)=='MOFA'){
  f_v<-get_factors(MOFAobject, factors = sel_factors[fn_sel] )[[1]]
  hist(f_v, breaks = 25)
  
  ws<-get_weights(MOFAobject, views = view, factors=sel_factors[fn_sel])[[1]]
  ws
  
  if (fn_sel==2){
    cut_high<-0.98; cut_low=0.02
    
  }else{
    cut_high<-0.998; cut_low=0.002
    
    
  } 
  cut_high<-0.5; cut_low=0.5
  
  ws_high<-ws[ws>quantile(ws, cut_high),]
  ws_low<-ws[ws<quantile(ws, cut_low),]
  ws_high
  
  ws_union<-c(ws_high, ws_low)
  #ws_high<-ws_high[1:15]
  length(ws_union)
  #if (view=='RNA'){
  # ens_genes<-rownames(assay(se_filt_V08))[grep(paste0(names(ws_union), collapse='|'), rownames(assay(se_filt_V08)))]
  #  feat_names=ens_genes
  
  #}else{
  feat_names=names(ws_union)
}


#}

feat_names_ens<-gsub('\\..*', '',feat_names)
feat_names_ens




## outliers

### choose from deseq@

#feat_names= sigLRT_genes$gene

#feat_names= sigLRT_genes$gene


merged_melt_orig<-merged_melt_orig_1[merged_melt_orig_1$variable %in% make.names(feat_names),]

levels(merged_melt_orig$variable)

ens<-gsub('\\..*', '',merged_melt_orig$variable)

if (view=='RNA'){
  symb<-get_symbols_vector(ens)
  merged_melt_orig$symbol<-symb
}else{
  merged_melt_orig$symbol<-merged_melt_orig$variable
  
}


#
### ### NOW match factors to samples
# CREATE GROUPS BY FACTOR 
############################################


sel_cohort=FALSE
sel_cohort<-c(1)



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
Z1_matched<-Z1[match(merged_melt$PATNO,patnos_z1) ]
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

groups_from_mofa_factors<-function(merged_melt, MOFAobject ){
  
  #'
  #' @param MOFAobject description
  #'

  ### cluster by one factor 
  groups_kmeans<-cluster_by_mofa_factors(MOFAobject, factors=sel_factors[fn_sel], centers=2)
  # OR CHOOSE ALL 
  groups_kmeans<-cluster_by_mofa_factors(MOFAobject, factors=sel_factors, centers=2)
  
#  Z1_matched<-Z1[match(merged_melt$PATNO,names(groups_kmeans$cluster)) ]
  
  

  high_label<-which.max(groups_kmeans$centers)
  #groups_kmeans$cluster<-ifelse(groups_kmeans$cluster==high_label, 'HighFactor', 'LowFactor')
  pats<-names(groups_kmeans$cluster)
  
  ### MATCH the clusters with the dataframe 
  kmeans_matched<-groups_kmeans$cluster[match(merged_melt$PATNO, pats )]
  kmeans_grouping<-factor(kmeans_matched)
  
  groups_from_mofa_factors
  
  
  return(kmeans_grouping)
  
}


### Important: create groups only for the patients.
merged_melt$kmeans_grouping<-groups_from_mofa_factors(merged_melt, MOFAobject)

merged_melt$grouping<-merged_melt$kmeans_grouping





na_ps<-unique(merged_melt[!is.na(merged_melt$kmeans_grouping),]$PATNO)
merged_melt_filt<-merged_melt[merged_melt$PATNO %in% na_ps, ]

### EDIT GROUP TO BE PATNO and GROUP!!
# OR JUST GROUP  
merged_melt_filt$VISIT



### Did they change scales???
### Breaks down into two groups based on grouping 
se_filt_V08_pd_g1<-se_filt_V08_pd[,se_filt_V08_pd$PATNO %in% levels(factor(merged_melt$grouping))[1] ]
se_filt_V08_pd_g2<-se_filt_V08_pd[,se_filt_V08_pd$PATNO %in% levels(factor(merged_melt$grouping))[2] ]


G1<-se_filt_V08_pd_g1[,!(se_filt_V08_pd_g1$PDSTATE == 'ON')]
G2<-se_filt_V08_pd_g2[,!(se_filt_V08_pd_g2$PDSTATE == 'ON')]
## TODO: plot here ALL the clinical variables by grouping!! 

#### 

merged_melt_filt$VISIT<-as.factor(merged_melt_filt$VISIT)

group_cats<-levels(factor(merged_melt$grouping))








################


### Plot to remove the other group ####
# TAKE THE low group  
# TODO: decide how to take the lowest x and highest x 
### TODO: DO THIS BOTH FOR CONTROLS AND DISEASE ####? 


merged_melt_filt$grouping<-merged_melt_filt$kmeans_grouping

merged_melt_filt$group<-as.logical(merged_melt_filt$grouping)
group_cat='grouping'
group_cat='Z2grouping'
group_cat='kmeans_grouping'

merged_melt_filt$group<-as.logical(merged_melt_filt[, group_cat])

merged_melt_filt$group<-as.factor(merged_melt_filt[, group_cat] )

merged_melt_filt_g1=merged_melt_filt[merged_melt_filt$group %in% group_cats[1],]
merged_melt_filt_g1=merged_melt_filt_g1[merged_melt_filt_g1$VISIT %in% c('BL', 'V08'),]
merged_melt_filt_g2=merged_melt_filt[merged_melt_filt$group %in% group_cats[2],]
merged_melt_filt_g2=merged_melt_filt_g2[merged_melt_filt_g2$VISIT %in% c('BL', 'V08'),]

merged_melt_filt_g1$VISIT<-as.factor(merged_melt_filt_g1$VISIT)
merged_melt_filt_g2$VISIT<-as.factor(merged_melt_filt_g2$VISIT)


######## First find out which of the molecules significantly change over time ####

#### TODO: do this ONLY  for disease AND SAVE THEM !! 
# THEREFORE MAKE THE GROUPING INTO A FUNCTION
# TODO: check both groups for significant changes
merged_melt_filt_g1


wilcox_stats1<-merged_melt_filt_g1 %>%
  group_by(symbol) %>%
  do(w=wilcox.test(value~VISIT, data=.)) %>%
  summarize(symbol, Wilcox=w$p.value) %>%
  as.data.frame()



wilcox_stats2<-merged_melt_filt_g2 %>%
  group_by(symbol) %>%
  do(w=wilcox.test(value~VISIT, data=.) )%>%
  summarize(symbol, Wilcox=w$p.value) %>%
  as.data.frame()


most_sig_over_time1<-wilcox_stats1[order(wilcox_stats1$Wilcox),][1:15,]
most_sig_over_time2<-wilcox_stats2[order(wilcox_stats2$Wilcox),][1:15,]

most_sig_over_time<-rbind(most_sig_over_time1, most_sig_over_time2)


#### CHOOSE 
merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in%  most_sig_over_time$symbol,]

merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in%  most_sig_over_time$symbol,]



### remove the ones insiude copntrols
wilcox_stats_controls<-merged_melt_filt_g1 %>%
  group_by(symbol) %>%
  do(w=wilcox.test(value~VISIT, data=.))%>%
  summarize(symbol, Wilcox=w$p.value) %>%
  as.data.frame()





most_sig_over_time1<-wilcox_stats1[order(wilcox_stats1$Wilcox),][1:15,]

most_sig_over_time2<-wilcox_stats2[order(wilcox_stats2$Wilcox),][1:15,]
#most_sig_over_time_controls<-wilcox_stats1_controls[order(wilcox_stats_controls$Wilcox),][1:15,]


# they should chgange in pd but not in controls!! 
# remove the ones in the controls
most_sig_over_time<-rbind(most_sig_over_time1, most_sig_over_time2)

most_sig_over_time_deseq = c('hsa.let.7a.3p', 'hsa.let.7f.1.3p', 'hsa.miR.101.3p', 'hsa.miR.142.5p')
####### CHOOSE 

#most_sig_over_time_deseq = make.names(sigLRT_genes$gene)
#most_sig_over_time_deseq<-make.names(colnames(data.filtered.only.pd))

#most_sig_over_time_deseq
merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in% most_sig_over_time_deseq[1:10],]


merged_melt_filt_g2_sig<-merged_melt_filt_g2[merged_melt_filt_g2$symbol %in%  most_sig_over_time$symbol,]
merged_melt_filt_g1_sig<-merged_melt_filt_g1[merged_melt_filt_g1$symbol %in%  most_sig_over_time$symbol,]




merged_melt_filt_g2_sig$COHORT=factor(merged_melt_filt_g2_sig$COHORT)
merged_melt_filt_g2_sig$VISIT=factor(merged_melt_filt_g2_sig$VISIT)



ggplot(data = merged_melt_all, aes(x = VISIT, y = value)) + 
  geom_point(aes(col=VISIT), size = 2) +
  geom_line(aes(group=PATNO),  col= 'grey') +
  geom_boxplot(aes(fill=VISIT))


### First answer : CAN THEY DIFFERENTIATE DISEASE CONTROL? 
# TODO: ADD DISEASE CONTROL

merged_melt_ct$kmeans_grouping='CONTROL'
merged_melt_ct$group='CONTROL'
merged_melt_ct$grouping='CONTROL'

filt_top=TRUE


### PUT THEM ALL TOGETHER IN THE BOXPLOTS 
#merged_melt_all<-rbind(merged_melt_ct, merged_melt_filt_g2_sig)
#merged_melt_all<-rbind(merged_melt_all, merged_melt_filt_g1_sig)



#### in the boxplots add the groups 
### first controls-- all markers need to be different in controls
### and second in the two groups of disease 
ggplot(data = merged_melt_filt_g2_sig, aes(x = VISIT, y = value, fill=kmeans_grouping)) + 
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
ggsave(paste0(outdir, '/trajectories/boxplots_', sel_factors[fn_sel],'_', view,'_',group_cat,sel_cohort , '.jpeg'), 
       width=12, height=12)


ggplot(data = merged_melt_filt_g1_sig, aes(x = VISIT, y = value, fill=COHORT)) + 
  #geom_point(aes(col=VISIT), size = 2) +
  #geom_line(aes(group=PATNO),  col= 'grey') +
  # subgroup should be in the fill parameter!!! 
  geom_boxplot(aes(x=VISIT, fill=COHORT ))+
  scale_color_viridis_d(option='mako')+
  scale_fill_viridis_d(option='mako')+
  
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +
  
  geom_signif(comparisons = list(c('BL', 'V08')),  
              map_signif_level=TRUE, 
              tip_length = 0, vjust=0.4)+
  
  facet_wrap(. ~ symbol, scales='free_y') +
  
  theme_bw() 
ggsave(paste0(outdir, '/trajectories/boxplots_', sel_factors[fn_sel],'_', view,'_',group_cat,sel_cohort , '.jpeg'), 
       width=12, height=12)




#######################################################
############ TIME TRAJECTORY FOR ALL VISITS ###########
#######################################################




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



colnames(merged_melt_filt)
colnames(merged_melt_ct)

merged_melt_ct$kmeans_grouping='CONTROL'
merged_melt_ct$group='CONTROL'
merged_melt_ct$grouping='CONTROL'


merged_melt_filt=rbind(merged_melt_filt,merged_melt_ct )

filt_top=TRUE


if (filt_top){
  
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% most_sig_over_time_deseq[1:10],]
  merged_melt_filt_most_sig<-merged_melt_filt[merged_melt_filt$symbol %in% most_sig_over_time$symbol[1:5],]
  
  
  nrow=NULL; height=2.6*4
}else{
  merged_melt_filt_most_sig<-merged_melt_filt
  nrow=NULL; height=7
  
}


merged_melt_filt_most_sig
### BY GROUP ####
#### TODO: plot also for CONTROLS! the same exact molecules thought.... so select them with PD 
## it only plots one group? 
ggplot(data = merged_melt_filt_most_sig, aes_string(x = 'VISIT', y = 'value', 
                                                    fill='group', group='group', colour='group')) + 
  stat_summary(geom = "pointrange", fun.data = median_IQR, 
               position=position_dodge(0))+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1) + 
  scale_color_viridis_d(option='turbo')+
  facet_wrap(. ~ symbol, scales='free_y', 
             nrow = nrow) +
  
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



#warnings()
ggsave(paste0(outdir, '/trajectories/trajectory', sel_factors[fn_sel],'_', view,  group_cat, filt_top,sel_cohort,  '.jpeg'), 
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
               geom = "line", size = 1) + 
### CHANGE OF MOLECULE VS CHANGE OF NP3



merged_melt_filt_most_sig$symbol
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
add_molecules_changes=FALSE
merged_melt_filt_most_sig$TOP=FALSE

if (add_molecules_changes){
  top_change<-molecules_change_by_patno[order(molecules_change_by_patno$diff, decreasing = FALSE)[1:20],'PATNO']
  top_change2<-molecules_change_by_patno[order(molecules_change_by_patno$log_FC, decreasing = FALSE)[1:20],'PATNO']
  just_molecules<-merged_melt[merged_melt$variable %in%  c(most_sig_over_time_deseq), ][, c('value', 'PATNO')]
  top_change3<-just_molecules[just_molecules$value< (-0),]$PATNO
  
  
  hist(merged_melt[merged_melt$variable %in%  c(most_sig_over_time_deseq), ]$value)
  top_change3
  merged_melt_filt_most_sig$TOP=FALSE
  merged_melt_filt_most_sig[merged_melt_filt_most_sig$PATNO %in% c( top_change2,top_change,top_change3) ,]$TOP<-TRUE
  any(merged_melt_filt_most_sig$TOP)
  
  
  top_molecular_patients<-c( top_change2,top_change,top_change3)
  top_molecular_patients
  
  
}




merged_melt_filt_most_sig



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
ggsave(paste0(outdir, '/trajectories/trajectory_by_pat_', sel_factors[fn_sel],'_', view,  group_cat, filt_top,'_', sel_cohort, '.jpeg'), 
       width=12, height=height)

merged_melt_filt_most_sig 

#df2$change

merged_melt$kmeans_grouping


#### ADD 18 month progression 

patnos_z1<-gsub('\\_.*', '', names(groups_kmeans$cluster))
groups_kmeans$cluster


Z2_grouping_df<-data.frame(group=groups_kmeans$cluster, PATNO=patnos_z1)
df18_months_2<-merge(df18_months, Z2_grouping_df, by='PATNO')
df18_months_2<-df18_months_2[!duplicated(df18_months_2),]
df18_months_2$value


ggplot(df18_months_2, aes(x=variable,y=value, group=group, colour=group)  )+
  geom_point(aes(x=variable,y=value, colour=group), alpha=0.5 )+
  # geom_line(aes(x=variable,y=value, group=PATNO, colour=group), lwd=0.2 )+
  stat_summary(fun = median, position=position_dodge(width=0), 
               geom = "line", size = 1.3) 




### CHANGE OF MOLECULE VS CHANGE OF NP3



merged_melt_filt_most_sig$symbol
levels(merged_melt_filt$symbol)
merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% 'hsa.miR.101.3p',]
merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% 'hsa.miR.101.3p',]
merged_melt_filt_1<-merged_melt_filt[merged_melt_filt$symbol %in% 'hsa.let.7a.3p',]


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

scale_change<-df_to_calc[,c( 'diff_scale', 'PATNO',paste0('PDSTATE_', sel_visit ))]
molecules_change_by_patno<-molecules_by_visit2[,c('log_FC','diff', 'PATNO', 'kmeans_grouping_V08')]

molecules_change_by_patno<-merge(molecules_change_by_patno, scale_change, by='PATNO')

molecules_change_by_patno[which.max(molecules_change_by_patno$diff),'PATNO']

hist(scale_change$diff_scale)


### color the top molecular ones too

scale_change$TOP=FALSE

scale_change[scale_change$PATNO%in%top_molecular_patients, ]$TOP=TRUE

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

ggplot(scale_change_gr, aes(x=diff_scale))+
  geom_histogram(aes(fill=TOP))


scale_change<-df_to_calc[,c( 'diff_scale', 'PATNO',paste0('PDSTATE_', sel_visit ))]


molecules_change_by_patno<-molecules_by_visit2[,c('log_FC','diff', 'PATNO', 'kmeans_grouping_V08')]
molecules_change_by_patno<-merge(molecules_change_by_patno, scale_change, by='PATNO')

# Plot the absolute difference between
# Diff
ggplot(molecules_change_by_patno, aes(x=diff, y=diff_scale))+
  geom_point(aes(color=kmeans_grouping_V08))+
  geom_smooth()+
  facet_wrap(~ PDSTATE_V16, nrow=3)


top_molecular_patients






##############

### Plot clinical variables ####

na_ps
### AND ALSO changed STAGE 

preprocess_visit_cl<-function(se_filt_V){
  
  se_filt_V_pd<-se_filt_V[,se_filt_V$COHORT == 1]
  se_filt_V_pd<-se_filt_V_pd[,se_filt_V_pd$PATNO %in% common]
  se_filt_V_pd<-se_filt_V_pd[,se_filt_V_pd$PATNO%in% na_ps ]
  
  return(se_filt_V_pd)
  
}
se_V08_cl<-preprocess_visit_cl(se_filt_V08)
se_V04_cl<-preprocess_visit_cl(se_filt_V04)
se_V06_cl<-preprocess_visit_cl(se_filt_V06)
se_bl_cl<-preprocess_visit_cl(se_filt_BL)



sel_factors[fn_sel]
cors_both<-get_correlations(MOFAobject, names(MOFAobject@samples_metadata))
cors_all=cors_both[[1]]
fn_sel
to_sel_cor<-names(cors_all[names(sel_factors)[fn_sel],][cors_all[names(sel_factors)[fn_sel],]>1.5])
to_sel<- c('PATNO', to_sel_cor)
to_sel
to_sel<-intersect(names(colData(se_V08_cl)), to_sel)
to_sel<-intersect(names(colData(se_bl_cl)), to_sel)

se_V08_cl$ESS_TOT
df_V08_cl<-colData(se_V08_cl)[,to_sel]
df_bl_cl<-colData(se_bl_cl)[,to_sel]
df_V06_cl<-colData(se_V06_cl)[,to_sel]
df_V04_cl<-colData(se_V04_cl)[,to_sel]

### Add also from other time points only the clinical data ####

df_V08_cl$PATNO






df_V08_cl_ml<-reshape2::melt(df_V08_cl,value.name = to_sel); df_V08_cl_ml$VISIT='V08'
df_bl_cl_ml<-reshape2::melt(df_bl_cl, value.name=to_sel); df_bl_cl_ml$VISIT='BL'
df_V04_cl_ml<-reshape2::melt(df_V04_cl, value.name=to_sel); df_V04_cl_ml$VISIT='V04'
df_V06_cl_ml<-reshape2::melt(df_V06_cl, value.name=to_sel); df_V06_cl_ml$VISIT='V06'




## Which variables are corelated? 


merged_melt_cl<-rbind(df_bl_cl_ml,df_V08_cl_ml)
merged_melt_cl<-rbind(merged_melt_cl,df_V06_cl_ml)
merged_melt_cl<-rbind(merged_melt_cl,df_V04_cl_ml)

dim(merged_melt_cl)

colnames(merged_melt_cl)



patnos_z1<-gsub('\\_.*', '', names(Z1))

Z1_matched<-Z1[match(merged_melt_cl$PATNO,patnos_z1) ]

hist(Z1_matched);
quantile(Z1_matched,0.9, na.rm=TRUE)
Z1_grouping<-factor(Z1_matched>quantile(Z1,0.9, na.rm=TRUE))
Z1_grouping<-factor(Z1_matched>quantile(Z1,0.8, na.rm=TRUE))

Z1_grouping<-factor(Z1_matched>quantile(Z1,0.8, na.rm=TRUE))






#####




sel_factors[fn_sel]
if (names(sel_factors[fn_sel]) %in% c('Factor4', 'Factor14', 'Factor1')){
  T=0.2
  Z1_grouping<-Z1_matched>quantile(Z1_matched,T, na.rm=TRUE)
  group_by_patient<-factor(Z1>quantile(Z1,0.2, na.rm=TRUE))
  
}

Z1_grouping
dim(merged_melt)

merged_melt_cl$grouping<-factor(ifelse(as.logical(Z1_grouping), 'HighFactor', 'LowFactor'))




### Create groups 

### Decide on the grouping #### 
### TODO: update for KMEANS grouping here too!! 
group_by_patient<-clusters_mofa$cluster


group_by_patient<-clusters$cluster
group_by_patient<-clusters_mofa_outcome$cluster
group_by_patient<-clusters_mofa$cluster
group_by_patient<-factor(Z1>quantile(Z1,0.75, na.rm=TRUE))


names(group_by_patient)<-gsub('\\_.*', '', names(group_by_patient))


combined_bl_log_common<-combined_bl_log[combined_bl_log$PATNO %in% na_ps,]

combined_bl_log_common

combined_bl_log_common$Z1_grouping<-group_by_patient[match(combined_bl_log_common$PATNO, names(group_by_patient))]



#combined_bl_log_common$grouping<-factor(ifelse(as.logical(combined_bl_log_common$Z1_grouping), 'HighFactor', 'LowFactor'))
combined_bl_log_common$grouping<-factor(ifelse(as.logical(combined_bl_log_common$Z1_grouping), 'HighFactor', 'LowFactor'))
combined_bl_log_common$grouping<-factor(combined_bl_log_common$Z1_grouping)
combined_bl_log_common$grouping

combined_bl_log_common$VISIT=factor(combined_bl_log_common$EVENT_ID)




### for categorical 
#table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3SPCH'],merged_melt_cl$grouping )
#table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3SPCH'], )
ggplot(data = combined_bl_log_common, aes( x=factor(grouping), 
                                           fill = factor(PDSTATE) )) + 
  geom_bar()+
  facet_wrap(. ~ VISIT, scales='free_y') 


ggplot(data = combined_bl_log_common, aes( x=factor(grouping), 
                                           fill = factor(PDSTATE) )) + 
  geom_bar()+
  facet_wrap(. ~ EVENT_ID, scales='free_y') 

#theme_bw() 
to_sel

### for continous 
is.numeric(merged_melt_cl$LAST_UPDATE_M1)


to_sel

merged_melt_cl$co
to_plot<-c('NP2PTOT','NP3TOT', 'NP3GAIT' , 'NP3BRADY', 'SCAU_TOT', 'scopa_cv', 
           'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP3_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT'
)

if (names(sel_factors[fn_sel]) %in% c('Factor3')){
  to_plot<-c('NP2PTOT','NP3TOT', 'NP3GAIT' , 'NP3BRADY', 'SCAU_TOT', 'scopa_cv', 
             'con_putamen', 'rigidity', 'td_pigd_old', 
             'RBD_TOT', 'NP3_TOT', 'NP2_TOT' , 'moca')
  
}else{
  to_plot<-c('NP2PTOT','NP3TOT' , 'NP3BRADY', 
             'td_pigd_old_on',  'AGE')
  to_plot<-c('NP2PTOT','NP3TOT', 'NP3GAIT' , 'NP3BRADY', 'SCAU_TOT', 'scopa_cv', 
             'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP3_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT' )
  to_plot<-selected_covars_broad
  
}

to_plot
## todo why is scau missing from baseline? how to measure total? 
merged_melt_cl_off<-merged_melt_cl[merged_melt_cl$PDSTATE %in% c('OFF', ''),]





### Filters #### 
# 1. 






remove_on=FALSE
if (remove_on){
  df_plot<- merged_melt_cl_off
  
  combined_bl_log_common_off=combined_bl_log_common[combined_bl_log_common$PDSTATE %in% c('OFF', ''),]
  
}else{
  df_plot<-merged_melt_cl
  
  
  combined_bl_log_common_off=combined_bl_log_common
  
}









combined_bl_log_common_off<-combined_bl_log_common_off[grep('BL|V04|V06|V08|V12|V16', combined_bl_log_common_off$EVENT_ID) ,]
#combined_bl_log_common_off<-combined_bl_log_common_off[grep('BL|V04|V06|V08', combined_bl_log_common_off$EVENT_ID) ,]

combined_bl_log_common_off<-combined_bl_log_common_off[grep('BL|V', combined_bl_log_common_off$EVENT_ID) ,]


df_plot<-combined_bl_log_common_off


df_plot %>% 
  group_by(EVENT_ID, grouping)%>% 
  summarize(count_distinct = n_distinct(PATNO))
time_nos<-df_plot %>% 
  group_by(EVENT_ID)%>% 
  summarize(count_distinct = n_distinct(PATNO))

smaller_group_id<-time_nos$EVENT_ID[which.min(time_nos$count_distinct)]
smaller_group<-df_plot[df_plot$EVENT_ID==smaller_group_id,]$PATNO

df_plot<-df_plot[df_plot$PATNO %in% smaller_group,]


df_V16<-df_plot[c(df_plot$VISIT=='V16'), ]
df_BL<-df_plot[df_plot$VISIT=='BL', ]

kruskal.test(df_V16$grouping,df_V16$NP3_TOT)
kruskal.test(df_BL$grouping,df_BL$NP3_TOT)


df_plot$grouping
######### boxplots by cluster over time
df_plot[,c('VISIT', 'grouping')][df_plot$VISIT=='V12',]
to_plot
df_plot_2k<-df_plot[df_plot$grouping %in% c(1,4),]

df_plot_2k<-df_plot[df_plot$grouping %in% c(4,5),]


df_plot_2k<-df_plot[df_plot$grouping %in% c(2,5),]
df_plot_2k<-df_plot

for (y in to_plot){
  
  ggplot(data = df_plot_2k, aes_string(x = 'VISIT', y = y, 
                                       fill='grouping', group='grouping', colour='grouping')) + 
    stat_summary(geom = "pointrange", fun.data = median_IQR, 
                 position=position_dodge(0))+
    stat_summary(fun = median, position=position_dodge(width=0), 
                 geom = "line", size = 1) + 
    scale_color_viridis_d(option='turbo')+
    #facet_wrap(. ~ symbol, scales='free_y', 
    #           nrow = 1) +
    
    #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
    theme_bw()+ 
    # geom_signif(comparisons = list(c('BL', 'V08')), 
    #            map_signif_level=TRUE, 
    #           tip_length = 0, vjust=0)+
    
    labs(y=y)+
    # legend(legend=c('Low', 'High'))+
    theme(strip.text = element_text(
      size = 10, color = "dark green"), 
      axis.title.y =element_text(
        size = 13, color = "dark green"), 
      axis.text.x = element_text(
        size = 9 ))
  
  
  
  warnings()
  ggsave(paste0(outdir, '/trajectories/trajectory_', sel_factors[fn_sel],'_', filt_top, y,'_', remove_on,  '.jpeg'), 
         width=5, height=3)
}


merged_melt_cl3<-merged_melt_cl
for (y in to_plot){
  ggplot(data = df_plot, aes_string(x = 'VISIT', y = y, 
                                    fill='grouping')) + 
    geom_boxplot()+
    scale_color_viridis_d(option='turbo')+
    #facet_wrap(. ~ symbol, scales='free_y', 
    #           nrow = 1) +
    
    #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
    theme_bw()+ 
    # geom_signif(comparisons = list(c('BL', 'V08')), 
    #            map_signif_level=TRUE, 
    #           tip_length = 0, vjust=0)+
    
    labs(y=y)+
    # legend(legend=c('Low', 'High'))+
    theme(strip.text = element_text(
      size = 10, color = "dark green"), 
      axis.title.y =element_text(
        size = 13, color = "dark green"), 
      axis.text.x = element_text(
        size = 9 ))
  
  
  
  warnings()
  ggsave(paste0(outdir, '/trajectories/box_', sel_factors[fn_sel],'_', filt_top, y,'_', remove_on,  '.jpeg'), 
         width=5, height=3)
}











merged_melt_cl3=merged_melt_cl
merged_melt_cl2=merged_melt_cl



for (cov_to_plot in to_sel){
  
  if (is.numeric(merged_melt_cl3[, cov_to_plot])){
    
    
    ggplot(data = merged_melt_cl3, aes_string(x = 'grouping', y =cov_to_plot  )) + 
      geom_boxplot(aes(fill=VISIT ))+
      #facet_wrap(. ~ symbol, scales='free_y') +
      theme_bw() 
    
    ggsave(paste0(outdir, '/trajectories/', 'cl_var', cov_to_plot, sel_factors[fn_sel],to_sel[2]  , '.jpeg'))
    merged_melt_cl2[,cov_to_plot]<-factor(merged_melt_cl2[,cov_to_plot])
    
    ggplot(data = merged_melt_cl2, aes_string( x='grouping', 
                                               fill = cov_to_plot) ) + 
      geom_bar()+
      facet_wrap(. ~ VISIT, scales='free_y') 
    
    ggsave(paste0(outdir, '/trajectories/', 'cl_var_discrete', cov_to_plot, sel_factors[fn_sel],to_sel[2]  , '.jpeg'))
    
  }
}

sm=samples_metadata(MOFAobject)


ggsave(paste0(outdir, '/trajectories/', sel_factors[fn_sel],to_sel[2]  , '.jpeg'))








table( merged_melt_cl[, 'VISIT'], merged_melt_cl[, 'NP3TOT'],merged_melt_cl$grouping )


ggplot(data = merged_melt_cl, aes_string(x = 'VISIT', y ='NP3SPCH'  )) + 
  geom_point(aes(col=factor(VISIT)), size = 2) +
  geom_line(aes(group=PATNO, col=VISIT), jitter()) +
  geom_jitter(aes(colour=grouping))+
  #  geom_boxplot(aes(fill=grouping))+
  #geom_line(aes(group=patno), palette='jco') +
  #facet_wrap(. ~ symbol) +
  
  #ggtitle(paste0('Factor ',sel_factors[fn_sel]))+
  theme_bw() 



#  stat_compare_means(comparisons = list(c("BL", "V08")), 
#                    label = "p.format", method = "wilcox.test", tip.length = 0)






### compare and group
library(ggpubr)
library(rstatix)









#dir.create(outdir, '/trajectories/')



#ggpaired(merged_melt, x = "VISIT", y = "value",
##        # color = "VISIT", line.color = "gray", line.size = 0.4,
#        color = "VISIT", line.color = 'grouping', line.size = 0.4,
#         
#         palette = "jco")+
# facet_wrap(facets = 'symbol', nrow=ceiling(length(ws_high)/3), ncol=4)+#
#facet_free(facets = 'symbol')+

#stat_compare_means(paired = TRUE)



merged_melt$patno


library(ggpubr)
# Create line plots of means
#ggline(merged_df, x = 'variable', y = "value", 
#       add = c("mean_sd", "jitter"))


