

##### ADD the groups from MOFA or other clusterings  #############
### Create groups 

### Decide on the grouping #### 
### TODO: update for KMEANS grouping here too!! 
# TODO: here take the grouping
# find all clustering types from mofa and sypply to function
group_by_patient<-clusters_mofa$cluster
group_by_patient<-clusters$cluster
group_by_patient<-clusters_mofa_outcome$cluster


########### HERE IT TAKES AS input all the metadata ############################
#names(group_by_patient)<-gsub('\\_.*', '', names(group_by_patient))


# TODO: HERE PLOT THE SCALES THAT ARE RELEVANT TO EACH FACTOR!!!!
## ir. check what are the corelations, and with which variables-for the PD patinets only

# where to get this from? 
all_fs_diff
all_diff_variables
#1. extract the diff variables and remove the diff




scale_vars_diff
imaging_variables_diff

to_plot<-c(scale_vars_diff)
## todo why is scau missing from baseline? how to measure total? 
# grep('BL|V04|V06|V08|V12|V16', combined_bl_log$EVENT_ID)
#either supply or grep 
### use df_mofa here?? since it is already set 


all_event_ids_p<-c('BL','V04','V06','V08','V10','V12','V14','V16', 'V18')
all_event_ids_p<-c('BL','V04','V06','V08','V10','V12','V14','V16', 'V18')

unique(curated_total[!is.na(curated_total$updrs2_score),'EVENT_ID'])



sm=MOFAobject@samples_metadata
### obtain all patient event ids to getr one row per patient!! 
patno_event_ids = sapply(all_event_ids_p, function(event_id){
                  return(paste0(sm$PATNO,'_', event_id ))
})
  
patno_event_ids=unlist(patno_event_ids)
# select data for the requested patiennts 
#combined_bl_log<-combined_new
combined_bl_log_sel<-fetch_metadata_by_patient_visit(patno_event_ids=patno_event_ids )
curated_mofa<-combined_bl_log_sel %>%
  dplyr::filter(EVENT_ID=='V14') %>%
 dplyr::filter(PATNO %in% sm$PATNO)


curated_mofa[, c('PATNO', 'updrs2_score')] %>%
      arrange(PATNO)

require(plyr)
combined_bl_log_sel$months <-as.numeric( mapvalues(combined_bl_log_sel$EVENT_ID, 
                               from= names(EVENT_MAP), 
                               to=unlist(EVENT_MAP, use.names=FALSE)))

lv_to_plot='V12'
combined_bl_log_sel$PDMEDYN_V14



########## Now we obtained the longitudinal input for these patients ####


# 1. Select the patients eg. by mofa groups 
df_to_attach<-combined_bl_log_sel
for (diff_var in names(all_clusts_mofa)){
  # 
  #   
        if (!is.null(all_clusts_mofa[[diff_var]])){
          
   
        clust_name = paste0(diff_var, '_clust')
         #print(clust_name)
         
        clusters_ids<-all_clusts_mofa[[diff_var]]
        df_to_attach[,clust_name]<-clusters_ids[match(df_to_attach$PATNO,names(clusters_ids ) )]
        
        df_to_attach[(df_to_attach$INEXPAGE %in% c('INEXHC')),paste0(diff_var, '_clust')]<-'HC'
        print(df_to_attach[,clust_name])
        }
}

combined_bl_log_sel<-df_to_attach

#combined_bl_log_common$grouping<-factor(ifelse(as.logical(combined_bl_log_common$Z1_grouping), 'HighFactor', 'LowFactor'))
combined_bl_log_sel$VISIT=factor(combined_bl_log_sel$EVENT_ID)
combined_bl_log_sel$month=factor(combined_bl_log_sel$months)

# TODO: Maybe get all time points? 

  
  
df_plot<-combined_bl_log_sel
## fetch grouping from MOFA 

df_plot_2k<-df_plot
df_plot_2k

PDSTATE_SEL=NULL

y='updrs2_score'
add_individual_lines=FALSE
add_boxplots<-TRUE
df_plot_2k<-df_plot
to_plot<-c('NP2PTOT','updrs2_score','updrs3_score','NP2PTOT', 'NP3TOT' ,'NP2_TOT','NP3_TOT', 'MCA_TOT', 'SCAU_TOT', 
           'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP1_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT' , 'scopa')




  
  ### TODO: check that only the variables with highest visit exist
  ###
  
  ## either loop through or melt 
  #df_plot_2k[df_plot_2k[,y]grouping]
  
  y='updrs2_score'
  
  df_plot_2k=df_plot
  

  
  plot_clinical_trajectory<-function(y, df_plot_2k=df_plot){
    
    
    lv='V13_V14';   
    lv='V12';  
    
    #clust_name<-paste0(y , '_diff_', lv,'_clust')
    #clust_name<-paste0(y , '_diff_', lv,'_clust')
    clust_name<-paste0('NP2PTOT' , '_clust')
    df_plot_2k<- df_plot_2k[df_plot_2k$months <=  EVENT_MAP[lv_to_plot],]
    clust_name %in% colnames(df_plot_2k)
    
    if (clust_name %in% colnames(df_plot_2k)){
      
      print(clust_name)
      
      df_plot_2k[, 'grouping']<-df_plot_2k[, clust_name]
      
      df_plot_2k$grouping<-as.factor(df_plot_2k$grouping)
      
      y_pl=y
      
      
      ### this is a function because for each subfilter other data is missing !! 
      df_plot_2k=df_plot_2k[!is.na(df_plot_2k$EVENT_ID),]
      df_plot_2k=df_plot_2k[!is.na(df_plot_2k[,clust_name]),]
      df_plot_2k[,y]=as.numeric( df_plot_2k[,y])
      
      df_lv<-df_plot_2k[df_plot_2k$EVENT_ID==lv_to_plot,]
      #df_lv<-df_plot_2k[df_plot_2k$EVENT_ID=='V12',]
      
      nums_plyr<-df_lv%>%
        dplyr::group_by(grouping) %>%
        dplyr::summarise(count = n_distinct(PATNO, grouping)) 
      

      #nums<-paste0('n=', paste0(table(unique(df_plot_2k[, c('grouping', 'PATNO')])[,'grouping'] ), collapse = ', '))
      nums<-paste0('n=', paste0(nums_plyr$count, collapse = ', '))
      

      p<-ggplot(data = df_plot_2k, aes_string(x = 'month', y = y, 
                                              fill='grouping',group='grouping',colour='grouping')) + 
        #stat_summary(geom = "pointrange", fun.data = median_IQR, 
        #             position=position_dodge(0), alpha=0.9)+
        stat_summary(fun = median, position=position_dodge(width=0), 
                     geom = "line", size = 1, alpha=0.7, lty='dashed', aes_string(colour='grouping')) 
      
      p
      #  geom_point(position='jitter', size=0.2)  
      # p
      #  stat_summary(fun = median_IQR, position=position_dodge(width=0), 
      #                   geom = "pointrange", size = 1, alpha=0.9) 
      
      if (add_individual_lines){
        p<-p+ geom_line(aes_string(x = 'month', y = y, 
                                   group='PATNO', colour='grouping' ),size=0.2, alpha=0.6)
      }
      
      if (add_boxplots){
        
        # p<-p+geom_violin(aes_string(x='VISIT', fill='grouping', group=NULL ), alpha=0.8)
        p<-p+geom_boxplot(aes_string(x='month', fill='grouping', group=NULL ),lwd=0.5, alpha=0.7)
        p
      }
      
      
      
      
      y_name=ifelse( (y %in% mt_kv[,1]), mt_kv[mt_kv[,1]==y,2], y)
      
      p<-p+scale_color_viridis_d(option='turbo')+
        scale_fill_viridis_d(option='turbo')+
        guides(fill=guide_legend(title='PD subgroup' ), color=guide_legend(title='PD subgroup' ))+
        
        
        labs(y=y_name, caption = paste0('group numbers: ', paste0(nums)))
      
      p<-p+ theme(axis.title.y =element_text(face='bold'))
      
      
      
      p
      warnings()
      ggsave(paste0(outdir, '/trajectories/clinical/trajectory_','_', y,'_',clust_name,'_lv_' , lv_to_plot, '.jpeg'), 
             width=5, height=3, dpi=300)
      
    }
    
    }
    
    
  # todo: are there duplicates in SOME patients only? does this change the medians and confidence intervals?
  for (y in to_plot){
    plot_clinical_trajectory(y, df_plot_2k=df_plot)
    
  }
 
  
    

  



#### ###################


## CLUSTER TRAJECTORIES ####

get_clinical_clusters(y)





combined_bl_log_sel
scale_mat=FALSE
y='NP3TOT'
y='NP3TOT'

y='SCAU_TOT'
y='stai_state'
y='stai_trait'
y='NP3TOT';
y='NP2PTOT'; nbCluster=3 #(last cluster is too small)
y='RBD_TOT'; nbCluster=3 #(last cluster is too small)


combined_bl_log_sel_pd$RBD_TOT==0

combined_bl_log_sel_pd$stai_state
    ### CLINICAcombined_bl_log_sel_pdL trajectories 
combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel$INEXPAGE=='INEXPD',]
combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel$INEXPAGE=='INEXPD',]

combined_bl_log_sel_pd<-combined_bl_log_sel_pd[!is.na(combined_bl_log_sel_pd$EVENT_ID),]
combined_bl_log_sel_pd[, y]=as.numeric(combined_bl_log_sel_pd[, y])
combined_bl_log_sel_pd<-combined_bl_log_sel_pd[!is.na(combined_bl_log_sel_pd$EVENT_ID),]


    df_plot_2k=combined_bl_log_sel_pd
    
    #y='NP2PTOT'
    
    x='EVENT_ID'
    df_plot_2k$PATNO=as.factor(df_plot_2k$PATNO)
    unique(combined_bl_log_sel_pd$PDMEDYN_V16)
    
    combined_bl_log_sel_pd_to_clust<-combined_bl_log_sel_pd[!(combined_bl_log_sel_pd$PDMEDYN_V14==0) ,]
    #combined_bl_log_sel_pd_to_clust<-combined_bl_log_sel_pd[combined_bl_log_sel_pd$PD_MED_USE_V12==1,]
    
    unique(combined_bl_log_sel_pd_to_clust$PATNO)
    cl_clusters_kml<-get_clinical_clusters_kml(combined_bl_log_sel_pd_to_clust,y, scale_mat = scale_mat, nbCluster=nbCluster )   
    cl_clusters<-cl_clusters_kml
    
    
    df_plot_2k$cluster<-as.factor(cl_clusters[match(df_plot_2k$PATNO,names(cl_clusters) )])
    df_plot_2k<-df_plot_2k[!is.na(df_plot_2k$cluster),]
    
    
    
    df_plot_2k_non_na<-df_plot_2k[!is.na(df_plot_2k[,y]), ]
    
    if (scale_mat){
      df_plot_2k_non_na[, y]<-scale(df_plot_2k_non_na[,y])
      
    }
    p<-ggplot(data = df_plot_2k_non_na, aes_string(x = x, y = y, colour='cluster', group='cluster'))
    p<-p+ 
      geom_line(aes_string(x = x, y = y ,group='PATNO', colour='cluster' ),size=0.2, alpha=0.3)+
      #geom_point(aes_string(x = x, y = y, colour='PDSTATE' ),size=0.6, alpha=0.6)+ 
      geom_point()+ 
      geom_smooth()
      #facet_wrap(~cluster)
    
    p
    
    
    ggsave(paste0(outdir, '/trajectories/clinical/clusters_',y, 'scale_',scale_mat, '.jpeg'), width=5, height=3)


    






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



