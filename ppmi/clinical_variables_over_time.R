

##### ADD the groups from MOFA or other clusterings  #############
### Create groups 

### Decide on the grouping #### 
### TODO: update for KMEANS grouping here too!! 
group_by_patient<-clusters_mofa$cluster


group_by_patient<-clusters$cluster
group_by_patient<-clusters_mofa_outcome$cluster

group_by_patient<- samples_metadata(MOFAobject)$N



########### HERE IT TAKES AS input all the metadata ############################
names(group_by_patient)<-gsub('\\_.*', '', names(group_by_patient))

# 1. Select the patients eg. by mofa groups 
combined_bl_log_common<-combined_bl_log[combined_bl_log$PATNO %in% na_ps,]
combined_bl_log_common$grouping<-group_by_patient[match(combined_bl_log_common$PATNO, names(group_by_patient))]



#combined_bl_log_common$grouping<-factor(ifelse(as.logical(combined_bl_log_common$Z1_grouping), 'HighFactor', 'LowFactor'))
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

# TODO: HERE PLOT THE SCALES THAT ARE RELEVANT TO EACH FACTOR!!!!
## ir. check what are the corelations, and with which variables-for the PD patinets only
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
  #to_plot<-selected_covars_broad
  ## note that these are calculated for future values!!! 
  to_plot<-c('NP2_TOT','NP3_TOT', 'MCA_TOT', 'SCAU_TOT', 
             'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP1_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT' )
  
}
to_plot<-c('NP2_TOT','NP3_TOT', 'MCA_TOT', 'SCAU_TOT', 
           'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP1_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT' )

merged_melt_cl$MCA_TOT
to_plot
## todo why is scau missing from baseline? how to measure total? 
merged_melt_cl_off<-merged_melt_cl[merged_melt_cl$PDSTATE %in% c('OFF', ''),]





### Filters #### 
# 1. 






PDSTATE_SEL='OFF'
df_plot<- merged_melt_cl_off
combined_bl_log_common_off=combined_bl_log_common[combined_bl_log_common$PDSTATE %in% c(PDSTATE_SEL),]

combined_bl_log_common_off$MCA_TOT

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
df_plot$grouping<-as.factor(df_plot$grouping)

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
    facet_wrap(~PDSTATE)+
    
    # legend(legend=c('Low', 'High'))+
    theme(strip.text = element_text(
      size = 10, color = "dark green"), 
      axis.title.y =element_text(
        size = 13, color = "dark green"), 
      axis.text.x = element_text(
        size = 9 ))
  
  
  
  warnings()
  ggsave(paste0(outdir, '/trajectories/clinical/trajectory_', factor,'_', filt_top, y,'_', PDSTATE_SEL,  '.jpeg'), 
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
    facet_wrap(~PDSTATE)+
    # legend(legend=c('Low', 'High'))+
    theme(strip.text = element_text(
      size = 10, color = "dark green"), 
      axis.title.y =element_text(
        size = 13, color = "dark green"), 
      axis.text.x = element_text(
        size = 9 ))
  
  
  
  warnings()
  ggsave(paste0(outdir, '/trajectories/clinical/box_', sel_factors[fn_sel],'_', filt_top, y,'_', PDSTATE_SEL,  '.jpeg'), 
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



