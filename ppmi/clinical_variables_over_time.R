

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
to_plot<-c('NP2_TOT','NP3_TOT', 'MCA_TOT', 'SCAU_TOT', 
           'con_putamen', 'rigidity', 'td_pigd_old', 'RBD_TOT', 'NP1_TOT', 'AGE_AT_VISIT', 'Outcome', 'NP4_TOT' )

to_plot<-c(scale_vars_diff)
## todo why is scau missing from baseline? how to measure total? 
# grep('BL|V04|V06|V08|V12|V16', combined_bl_log$EVENT_ID)
#either supply or grep 
### use df_mofa here?? since it is already set 
all_event_ids<-c('BL','V02','V04','V06','V08','V10','V12','V14','V16', 'V18')
all_event_ids<-unique(combined_bl_log$EVENT_ID)
all_event_ids_p<-c('BL',all_event_ids[grep('V', all_event_ids)] )
all_event_ids_p<-c('BL','V02','V04','V06','V08','V10','V12','V14','V16', 'V18')

sm=MOFAobject@samples_metadata
### obtain all patient event ids to getr one row per patient!! 
patno_event_ids = sapply(all_event_ids_p, function(event_id){
                  return(paste0(sm$PATNO,'_', event_id ))
})
  
patno_event_ids=unlist(patno_event_ids)
# select data for the requested patiennts 
#combined_bl_log<-combined_new
combined_bl_log_sel<-fetch_metadata_by_patient_visit(patno_event_ids=patno_event_ids )
unique(combined_bl_log_sel$EVENT_ID)
combined_bl_log_sel$PDMEDYN_V14
########## Now we obtained the longitudinal input for these patients ####


# 1. Select the patients eg. by mofa groups 
clust_y=all_clusts[1];
clust_y_labs_all<-sapply(all_clusts, function(clust_y){
      if (clust_y %in% colnames(MOFAobjectPD@samples_metadata)){
        group_by_patient<- MOFAobjectPD@samples_metadata[, clust_y]
        names(group_by_patient)<- MOFAobjectPD@samples_metadata$PATNO
        clust_y_labs<-group_by_patient[match(combined_bl_log_sel$PATNO, names(group_by_patient))];
        
        #combined_bl_log_sel[, clust_y]<-clust_y_labs
      print('attached')
      return(clust_y_labs)
        }
      
})
combined_bl_log_sel<-cbind(combined_bl_log_sel,clust_y_labs_all );

#combined_bl_log_common$grouping<-factor(ifelse(as.logical(combined_bl_log_common$Z1_grouping), 'HighFactor', 'LowFactor'))
combined_bl_log_sel$VISIT=factor(combined_bl_log_sel$EVENT_ID)
# TODO: Maybe get all time points? 

  
  
df_plot<-combined_bl_log_sel
## fetch grouping from MOFA 

df_plot_2k<-df_plot
df_plot_2k

PDSTATE_SEL=NULL

y='NP2PTOT'
add_individual_lines=FALSE
add_boxplots<-FALSE

for (y in to_plot){
  
  ### TODO: check that only the variables with highest visit exist
  ###
  
  ## either loop through or melt 
  #df_plot_2k[df_plot_2k[,y]grouping]
  lv='V16'
  clust_name<-paste0(y , '_diff_', lv,'_clust')
  if (clust_name %in% colnames(df_plot_2k)){
    df_plot_2k[, 'grouping']<-df_plot_2k[, clust_name]
    df_plot_2k$grouping<-as.factor(df_plot_2k$grouping)
    
  y_pl=y

  

  df_plot_2k=df_plot_2k[!is.na(df_plot_2k$EVENT_ID),]
  df_plot_2k=df_plot_2k[!is.na(df_plot_2k[,clust_name]),]
  
  df_lv<-df_plot_2k[df_plot_2k$EVENT_ID=='V16',]
  nums<-df_lv%>%
    group_by(grouping) %>%
    summarise(count = n_distinct(PATNO)) 
  
  
  p<-ggplot(data = df_plot_2k, aes_string(x = 'VISIT', y = y, 
                                       fill='grouping',group='grouping',  colour='grouping')) + 
    stat_summary(geom = "pointrange", fun.data = median_IQR, 
                 position=position_dodge(0), alpha=0.9)
    
    if (add_individual_lines){
     p<-p+ geom_line(aes_string(x = 'VISIT', y = y, 
                           group='PATNO', colour='grouping' ),size=0.2, alpha=0.6)
    }
  
    if (add_boxplots){

        p<-p+geom_violin(aes_string(x='VISIT', fill='grouping', group=NULL ), alpha=0.8)
        
    }
    
    
    
    p<-p+stat_summary(fun = median, position=position_dodge(width=0), 
                 geom = "line", size = 1, alpha=0.9) + 
    scale_color_viridis_d(option='magma')+
 
    # geom_signif(comparisons = list(c('BL', 'V08')), 
    #            map_signif_level=TRUE, 
    #           tip_length = 0, vjust=0)+
    
    labs(y=y, caption = paste0('group numbers: ', paste0(nums$count, collapse=', ')))

    theme(strip.text = element_text(
      size = 10, color = "dark green"), 
      axis.title.y =element_text(
        size = 13, color = "dark green"), 
      axis.text.x = element_text(
        size = 9 ))
  
  
  p
  warnings()
  ggsave(paste0(outdir, '/trajectories/clinical/trajectory_', factor,'_', filt_top, y,'_',clust_name,  '.jpeg'), 
         width=5, height=3)
  
  }
}


#### ###################


## CLUSTER TRAJECTORIES ####


get_clinical_clusters<-function(y, centers=4){
  #'
  #' @param 
  #'
  #'
        combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel$COHORT==1,]
        clin_traj<-combined_bl_log_sel[,c('PATNO','EVENT_ID', y)]
        
        clin_traj<-clin_traj[!is.na(clin_traj$EVENT_ID),]
        #clin_traj$months<-unlist(EVENT_MAP[clin_traj$EVENT_ID], use.names = FALSE)
        
        clin_traj_wide<-reshape(clin_traj, idvar='PATNO', timevar='EVENT_ID', direction='wide')
        rownames(clin_traj_wide)<-clin_traj_wide$PATNO
        clinical_clusters<-kmeans((na.omit(clin_traj_wide)[, -1]), centers=centers)
        
        return(clinical_clusters$cluster)
  
}
get_clinical_clusters(y)



get_clinical_clusters_kml<-function(combined_bl_log_sel_pd,y, nbCluster=4, scale_mat=FALSE){
  
  #combined_bl_log_sel_pd=combined_bl_log_sel_pd_to_clust
  #combined_bl_log_sel_pd<-combined_bl_log_sel[combined_bl_log_sel[,'INEXPAGE']=='INEXPD',]
  unique(combined_bl_log_sel_pd$EVENT_ID)
  
  clin_traj<-combined_bl_log_sel_pd[,c('PATNO','EVENT_ID', y)]
  
  clin_traj<-clin_traj[!is.na(clin_traj$EVENT_ID),]
  
  unique(clin_traj$EVENT_ID)
  
  
  clin_traj_wide<-reshape(clin_traj, idvar='PATNO', timevar='EVENT_ID', direction='wide')
  rownames(clin_traj_wide)<-clin_traj_wide$PATNO
  #clinical_clusters<-kmeans((na.omit(clin_traj_wide)[, -1]), centers=centers)
  
  #return(clinical_clusters$cluster)
  
  #install.packages('kml')
  
  ### Clinical trajectory means 
  # REMOVE columns full of NA
  clin_traj_mat<-as.data.frame(sapply((clin_traj_wide)[, -1], as.numeric))
  # ALSO SCALE
  #clin_traj_mat<-clin_traj_mat[!is.na(clin_traj_mat$EVENT_ID),]
  df<-clin_traj_mat
  clin_traj_mat <- as.matrix(df[,colSums(is.na(df))<nrow(df)])
  if (scale_mat){
    clin_traj_mat<-scale(clin_traj_mat)
    
  }
  #devtools::install_github("JimMcL/trajr")
  #library('trajr')
  #trj <- TrajGenerate(200, random = TRUE, angularErrorSd = .25)
  
  #smoothed<-TrajSmoothSG(clin_traj_mat[1,],3,31)
  
  
  
  CLD <- kml::cld(clin_traj_mat, timeInData = 1:dim(clin_traj_mat)[2], maxNA = 2)
  #length(CLD)
  clusters<-kml::kml(CLD, nbRedrawing = 5)
  
  #nbCluster=4
  
  # run choice
  clust_ids<-getClusters(CLD,nbCluster=nbCluster )
  names(clust_ids)<-clin_traj_wide$PATNO
  
  length(clust_ids)
  return(clust_ids)
}


library('kml')
library('dplyr')

combined_bl_log_sel
scale_mat=FALSE
y='NP3TOT'
y='NP3TOT'

y='SCAU_TOT'
y='stai_state'
y='stai_trait'
y='NP3TOT';
y='NP2PTOT'; nbCluster=3 #(last cluster is too small)


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



