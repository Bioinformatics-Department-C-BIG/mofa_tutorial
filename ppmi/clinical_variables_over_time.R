
##### ADD the groups from MOFA or other clusterings  #############
#1. extract the diff variables and remove the diff
to_plot<-c(scale_vars_diff)
all_event_ids_p<-c('BL','V04','V06','V08','V10','V12','V14','V16', 'V18')
sm=MOFAobject@samples_metadata
### obtain all patient event ids to get one row per patient!! 
patno_event_ids = sapply(all_event_ids_p, function(event_id){
                  return(paste0(sm$PATNO,'_', event_id ))
})
  
patno_event_ids=unlist(patno_event_ids)
# select data for the requested patiennts 
#combined_bl_log<-combined_new
combined_bl_log_sel<-fetch_metadata_by_patient_visit(patno_event_ids=patno_event_ids ) # todo - filter selection on or off 
curated_mofa<-combined_bl_log_sel %>%
  dplyr::filter(EVENT_ID=='V14') %>%
 dplyr::filter(PATNO %in% sm$PATNO)

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

df_plot<-combined_bl_log_sel
## fetch grouping from MOFA 

df_plot_2k<-df_plot


df_plot
PDSTATE_SEL=NULL

y='updrs2_score'
add_individual_lines=FALSE
add_boxplots<-FALSE
to_plot<-c('NP2PTOT','updrs2_score','updrs3_score','NP2PTOT', 'NP3TOT' ,'NP2_TOT','NP3_TOT', 'MCA_TOT', 'SCAU_TOT', 
           'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP1_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT' , 'scopa')
y='updrs2_score'
  
  
lv='V13_V14';   

  
  plot_clinical_trajectory<-function(y, clust_name, df_plot_2k=df_plot, lv='V12', fname=fname, add_individual_lines=FALSE){
    #'
    #' Plot the clinical value trajectory over time, until the last visit defined for each subgroup
    #' takes the clinival score, the last visit, and the cluster ids 
    #' @param y clinical value to plot
    #' @param y
    #' TODO: adjust the function to take grouping as a variable either from molecular or from clinical clusters
    #'  Created parameters: clust_name, grouping variable,


    df_plot_2k<- df_plot_2k[df_plot_2k$months <=  EVENT_MAP[lv_to_plot],]
    
    if (clust_name %in% colnames(df_plot_2k)){
      print(clust_name)
      df_plot_2k[, 'grouping']<-as.factor(df_plot_2k[, clust_name])
      y_pl=y
      ### filter inside the function because for each subfilter other data is missing !! 
      df_plot_2k=df_plot_2k[!is.na(df_plot_2k$EVENT_ID),]
      df_plot_2k=df_plot_2k[!is.na(df_plot_2k[,clust_name]),]
      df_plot_2k[,y]=as.numeric( df_plot_2k[,y])
      df_lv<-df_plot_2k[df_plot_2k$EVENT_ID==lv_to_plot,]
    
      # get frequencies in each cluster 
      # TODO: create function...? this is used also in boxplots and other metric calculations
      # it might be differnt per metric though so we need to recalculate
      nums_plyr<-df_lv%>%
        dplyr::group_by(grouping) %>%
        dplyr::summarise(count = n_distinct(PATNO, grouping)) 
      nums<-paste0('n=', paste0(nums_plyr$count, collapse = ', '))
    


      p<-ggplot(data = df_plot_2k, aes_string(x = 'month', y = y,
                fill='grouping',group='grouping',colour='grouping')) + 
        stat_summary(geom = "errorbar", fun.data = median_IQR, 
                     position=position_dodge(0.8), alpha=0.9)+
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
      ggsave(paste0(fname,  '.jpeg'), 
             width=5, height=3, dpi=300)
      
    }
    
    }
    
add_individual_lines=TRUE
to_plot=c('updrs2_score', 'NP2PTOT', 'MCATOT', 'moca')
clust_metric='moca'
clust_metric='NP2PTPT_LOG'

all_fs_diff$updrs2_score_LOG
df_plot
  # todo: are there duplicates in SOME patients only? does this change the medians and confidence intervals?
  for (y in to_plot){

    clust_name<-paste0(clust_metric, '_clust')

    cluster_params<-paste0(clust_name ,'/', k_centers_m,'/',rescale_option)
    cluster_params
    fname=paste0(outdir, cluster_params, '/trajectories/clinical/traj_','_', y,'_','_lv_' , lv_to_plot)
    plot_clinical_trajectory(y,clust_name=clust_name, df_plot_2k=df_plot,  lv='V12', fname=fname, add_individual_lines=add_individual_lines)
  }

  
    

  



#### ###################


## CLUSTER TRAJECTORIES ####
y='NP2PTOT'
get_clinical_clusters(y)


