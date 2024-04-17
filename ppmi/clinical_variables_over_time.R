
##### ADD the groups from MOFA or other clusterings  #############
#1. extract the diff variables and remove the diff
to_plot<-c(scale_vars_diff)
all_event_ids_p<-c('BL','V04','V06','V08','V10','V12','V14','V16', 'V18')
all_event_ids_p<-c('BL','V08','V10','V12','V14','V16', 'V18')

all_event_ids_match_molecular<-c('BL','V04','V06' ,'V08')


## metadata to select PATNOs FROM 

mofa_filter<-TRUE
if (mofa_filter){
  sm=MOFAobject@samples_metadata
}else{
  sm=combined_bl_log
}


sm_sel=sm



### obtain all patient event ids to get one row per patient! 
patno_event_ids = unlist(sapply(all_event_ids_p, function(event_id){
                  return(paste0(sm_sel$PATNO,'_', event_id ))
}))

patno_event_ids_mol = unlist(sapply(all_event_ids_match_molecular, function(event_id){
                  return(paste0(sm_sel$PATNO,'_', event_id ))
}))

# select data for the requested patiennts 
#combined_bl_log<-combined_new
combined_bl_log_sel<-fetch_metadata_by_patient_visit(patno_event_ids=patno_event_ids ) # todo - filter selection on or off 

combined_bl_log_sel_mol<-fetch_metadata_by_patient_visit(patno_event_ids=patno_event_ids_mol ) # todo - filter selection on or off 



curated_mofa<-combined_bl_log_sel %>%
  dplyr::filter(EVENT_ID=='V14') %>%
 dplyr::filter(PATNO %in% sm$PATNO)

require(plyr)
combined_bl_log_sel$months <-as.numeric( mapvalues(combined_bl_log_sel$EVENT_ID, 
                               from= names(EVENT_MAP_YEAR_NUM), 
                               to=unlist(EVENT_MAP_YEAR_NUM, use.names=FALSE)))

lv_to_plot='V12'
combined_bl_log_sel$PDMEDYN_V14

estimations_matched_all_combined<-estimations[match(combined_bl_log_sel$PATNO_EVENT_ID, rownames(estimations) ),]
combined_bl_log_sel<-cbind(combined_bl_log_sel,estimations_matched_all_combined)


########## Now we obtained the longitudinal input for these patients ####


# 1. Select the patients eg. by mofa groups 
df_to_attach<-combined_bl_log_sel

names(all_clusts_mofa)

df_to_attach<-attach_cluster_ids(df_to_attach, all_clusts_mofa)
combined_bl_log_sel_mol<-attach_cluster_ids(combined_bl_log_sel_mol, all_clusts_mofa)
combined_bl_log_sel_mol$NP3TOT_clust

combined_bl_log_sel<-df_to_attach




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

  lv_to_plot = 'V08'
  clust_name='COHORT_DEFINITION'
    clust_name='INEXPAGE'
  

times_sel = c('BL','V04', 'V06', 'V08')

cluster='NP3TOT_clust'

df_plot_mol<-combined_bl_log_sel_mol

combined_bl_log_sel_mol$NP3TOT


clust_name
medians_all_clusts<-get_variables_by_cluster_all_time(combined_bl_log_sel_mol, paste0(DIFF_VAR, '_clust'))
medians_all_clusts$cluster

medians_all_clusts

     # rownames(means_by_cluster)<-means_by_cluster$cluster









  plot_clinical_trajectory<-function(y, clust_name, df_plot_2k=df_plot, lv='V12', fname=fname, 
  add_individual_lines=FALSE, add_boxplots=FALSE, pal='turbo'){
    #'
    #' Plot the clinical value trajectory over time, until the last visit defined for each subgroup
    #' takes the clinival score, the last visit, and the cluster ids 
    #' @param y clinical value to plot
    #' @param y
    #' TODO: adjust the function to take grouping as a variable either from molecular or from clinical clusters
    #'  Created parameters: clust_name, grouping variable,
    #' 

    df_plot_2k<- df_plot_2k[df_plot_2k$months <=  EVENT_MAP_YEAR_NUM[lv_to_plot],]
    df_plot_2k$PATNO %in% samples_metadata(MOFAobject)$PATNO

    
    if (clust_name %in% colnames(df_plot_2k)){
      #print(clust_name)


      df_plot_2k[, 'grouping']<-as.factor(df_plot_2k[, clust_name])
      y_pl=y
      ### filter inside the function because for each subfilter other data is missing !! 
      df_plot_2k=df_plot_2k[!is.na(df_plot_2k$EVENT_ID),]
      df_plot_2k=df_plot_2k[!is.na(df_plot_2k[,clust_name]),]
      df_plot_2k=df_plot_2k[!df_plot_2k[,clust_name]=='',] # not empty
      

      table(df_plot_2k[,clust_name])

      df_plot_2k[,y]=as.numeric( df_plot_2k[,y])
      df_lv<-df_plot_2k[df_plot_2k$EVENT_ID==lv_to_plot,]
    
      # get frequencies in each cluster 
      # TODO: create function...? this is used also in boxplots and other metric calculations
      # it might be differnt per metric though so we need to recalculate
      nums_plyr<-df_lv%>%
        dplyr::group_by(grouping) %>%
        dplyr::summarise(count = n_distinct(PATNO, grouping)) 
      nums<-paste0('n=', paste0(nums_plyr$count, collapse = ', '))
    

#df_plot_2k$B.Memory
   # y='B.Memory'
      p<-ggplot(data = df_plot_2k, aes_string(x = 'month', y = y,
                fill='grouping',group='grouping',colour='grouping'))+

                
                 stat_summary(fun = median, position=position_dodge(width=0), 
                     geom = "line", size = 1, alpha=0.7, lty='dashed', aes_string(colour='grouping')) 

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
        p<-p+geom_boxplot(aes_string(x='month', fill='grouping', group=NULL ),lwd=0.2, alpha=0.7, 
        outlier.shape = NA)
        #scale_y_continuous(expand = expansion(mult = c(0.05, 0.005)))
        p<-p+ geom_pwc( tip.length = 0,
            method = "wilcox_test", label = "p.adj.signif", label.size = 2,
              bracket.nudge.y = -0.3)
        p
      }else{
         p<-p+
        stat_summary(geom = "errorbar", fun.data = median_IQR, # only add errorbar ifno box plot 
                    position=position_dodge(0.8), alpha=0.9)
       
      
      p
      }

      y_name=ifelse( (y %in% mt_kv[,1]), mt_kv[mt_kv[,1]==y,2], y)
      
      p<-p+scale_color_viridis_d(option=pal)+
        scale_fill_viridis_d(option=pal)+
        guides(fill=guide_legend(title='PD subgroup' ), color=guide_legend(title='PD subgroup' ))+
        
        
        labs(y=y_name,x='year', caption = paste0('group numbers: ', paste0(nums)))
      
      p<-p+ theme(axis.title.y =element_text(face='bold'))

      p
      warnings()
      


      width=1+log(1+as.numeric(gsub('V','',lv_to_plot)))*log(length(unique(df_plot_2k$grouping)))
      height=3

      ggsave(paste0(fname,  '.jpeg'), 
             width=width, height=height, dpi=300)
      
    }
    
    }

dir.create(paste0(cluster_params_dir, '/tr/clinical/'), recursive = TRUE)


## Plot cell types by COHORT 
to_plot<-c(colnames(estimations), 'Lymphocytes....', 'Neutrophils....')
for (y in  to_plot){
  df_plot_2k =df_plot_2k[!df_plot_2k$COHORT==4,]
  df_plot_2k =df_plot_2k[df_plot_2k$INEXPAGE %in% c(sel_subcoh, 'INEXHC'),]

  fname=paste0(outdir,'/tr/', clust_name, '_',  y, lv_to_plot  )
  plot_clinical_trajectory(y, clust_name, df_plot_2k, lv='V12', fname, add_boxplots = TRUE, pal='turbo' )

}

df_plot_2k$NP2PTOT_LOG_clust



## Plot cell types by Cluster! 
add_individual_lines=FALSE
sm$T.CD4.Naive
to_plot=c('updrs2_score', 'NP2PTOT', 'MCATOT', 'moca')
to_plot=c( 'Lymphocytes....', 'Neutrophils....', 'T.CD4.Naive', 'B.memory')

lv_to_plot = 'V08'
to_plot=c( 'Lymphocytes....', 'Neutrophils....','Neutrophils.LD', 'T.CD4.Naive', 'B.Memory', 'T.CD8.Memory')
lv_to_plot = 'V14'

to_plot=c('NP3TOT','updrs2_score','updrs3_score_on', 'NP2PTOT', 'MCATOT', 'moca', 'abeta', 'sft')
clust_metric='updrs3_score_on'
clust_metrics=c('moca', 'NP2PTOT_LOG', 'NP3TOT_LOG', 'updrs3_score_on')


  # todo: are there duplicates in SOME patients only? does this change the medians and confidence intervals?
  for (clust_metric in clust_metrics){
    for (y in to_plot){

      clust_name<-paste0(clust_metric, '_clust')
      fact<-get_factors_for_metric(clust_metric); fact_s=paste0(fact, collapse='_')
      cluster_params<-paste0( '/clustering/',  fact_s,'/', k_centers_m,'/r',as.numeric(rescale_option), '/')


      dir.create(paste0(outdir, cluster_params, '/tr/clinical/'), recursive = TRUE)

      fname=paste0(outdir, cluster_params, '/tr/clinical/traj_','_', y,'_','_lv_' , lv_to_plot)
      plot_clinical_trajectory(y,clust_name=clust_name, df_plot_2k=df_plot,  lv='V08', fname=fname,
      add_individual_lines=add_individual_lines, 
      add_boxplots = TRUE)
    
    
    }

    graphics.off()
  }

  
fname


#### ###################


## CLUSTER TRAJECTORIES ####
y='NP2PTOT'
get_clinical_clusters(y)


























