










sel_visit='V16'




get_future_clinvars<-function(combined_bl_log){
  
  
  combined_by_visit<-split(combined_bl_log, combined_bl_log$EVENT_ID )
  # 72 MONTHS 
  cm2<-combined_by_visit[c('BL','V08','V12', 'V14'  ,'V16','V18')]
  cm2<-combined_by_visit[c('BL','V08','V14')]
  cm2<-combined_by_visit[c('BL','V08','V12', 'V14'  ,'V16')]
  
  # cm2<-combined_by_visit[c('BL' ,'V12')]
  
  combinded_wide_all<-cm2 %>% 
    imap(function(x, y) x %>% rename_with(~paste(., y, sep = '_'), -PATNO)) %>%
    reduce(full_join, by = "PATNO")
  
  
  
  combinded_wide_all$np3tot
  combinded_wide_all[combinded_wide_all$COHORT_V16==1,]$scopa_BL
  combinded_wide<-combinded_wide_all[!is.na(combinded_wide_all[, paste0( cl_var,'_', sel_visit)]),]
  df2<-combinded_wide
  
  # df2<-df2%>% dplyr::filter(PDSTATE.y%in%c('OFF', ' '))
  #df2<-df2[,!colnames(df2)=='PDSTATE.y']
  df2<-df2[!duplicated(df2),]
  
  
  id_vars<-colnames(df2)[!grepl( cl_var, colnames(df2))]
  # id_vars<-colnames(df2)[!grepl( 'NP3TOT|PDSTATE', colnames(df2))]
  
  df2_melt<-reshape::melt(df2, id=c( id_vars))
  df2$COHORT_BL
  
  colnames(df2_melt)
  #df2_melt=df2_melt[df2_melt$COHORT_BL==1 & df2_melt$PATNO %in% sel_pats,]
  # TODO: should we filter with specific patients ? 
  sel_pats<-MOFAobject@samples_metadata$PATNO
  #df2_melt=df2_melt[ df2_melt$PATNO %in% sel_pats,]
  df2_melt=df2_melt[ df2_melt$PATNO %in% sel_pats,]
  
  
  #df2_melt$variable=as.numeric(df2_melt$variable)
  df_future_clinvars<-df2_melt

  
  
  ### Clinical scale trends per patient - is there a difference with controls ?
  # ISSUE WE CANNOT SEE IT for controls -it was not measured..
  graphics.off()
  pp<-ggplot(df2_melt, aes(x=variable,y=value)  )+
    geom_point(aes(x=variable,y=value, color=PATNO), size=0.2 )+
    geom_line(aes(x=variable,y=value, group=PATNO, color=PATNO), lwd=0.3, alpha=0.4) +
    scale_color_viridis_c(option='turbo')+
    facet_wrap(~COHORT_BL, nrow=2)
  
  return(df_future_clinvars)
}


clip_outliers_times<-function(X, x_times=1.5, upper=TRUE, lower=TRUE){
  #'
  #'
  #' outliers more than x-mes the IQR to be cliped 
  #'
  #'
  quartiles <- quantile(X, probs=c(.25, .75), na.rm = TRUE)
  IQR <- IQR(X, na.rm=TRUE)
  # Q3 + 1.5*IQR
  Lower <- quartiles[1] - x_times*IQR
  Upper <- quartiles[2] + x_times*IQR 
  
  if (upper){
    X[X>Upper]<-Upper
    
  }
  if (lower){
    
    X[X<Lower]<-Lower
  }
  return(X)
  
}




get_clinvar_changes<-function(df_future_clinvars, sel_visit, cl_var, sel_state){
  # TODO: combine multiple clinical variables 
  
  #' Calculate change
  #' @param
  #'
  #'
  
  
  
  #  facet_wrap('', nrow=3)
  
  ## filter for off:
  df2_off<-df_future_clinvars
  df2_off<-df2[df2[, paste0('PDSTATE_',sel_visit)]==sel_state, ]
  

  df_to_calc<-df2_off
  
  
  ### HOW many are controls? 
  #table(df2[df2$PATNO %in% sel_pats,]$COHORT_BL)
  #table(unique(df2_melt[,c('PATNO', 'PDMEDYN_V14', 'PDSTATE_V14', 'DBSYN_V14') ])[c('PDMEDYN_V14','DBSYN_V14' )])
  #table(unique(df2_melt[,c('PATNO', 'PDSTATE_V14') ])$PDSTATE_V14)
  #table(unique(df2_melt[,c('PATNO', 'PDSTATE_V14') ][,c('PDMEDYN_V14','DBSYN_V14' )]))
  
  df_to_calc$log_FC_scale<-log2(log2(df_to_calc[,paste0(cl_var,'_',sel_visit)])/ log2(df_to_calc[, paste0(cl_var, '_BL')]))
  
  X2_cl=df_to_calc[,paste0(cl_var,'_',sel_visit)]
  X1_cl=df_to_calc[,paste0(cl_var,'_','BL')]
  
 
  X2_cl<-clip_outliers(X2_cl)
  
  
  calc_change<-function(X1,X2){
    change<-((X2-X1)/(X2+X1))
  }
  calc_change2<-function(X1,X2){
    change<-log2(X2/X1)
  }
  calc_change_diff<-function(X1,X2){
    change<-X2-X1
  }
  
  
  df_to_calc$FC_scale<-calc_change(X1_cl, X2_cl)
  df_to_calc$log_FC_scale<-calc_change2(X1_cl, X2_cl)
  df_to_calc$diff_scale<-calc_change_diff(X1_cl, X2_cl)
  
  #df_to_calc$log_FC<-log2(log(X2)/log(X1))
  df_to_calc$COHORT_BL=factor(df_to_calc$COHORT_BL)
  
  median(df_to_calc$log_FC, na.rm=TRUE)
  return(df_to_calc)
}





######
#merged_melt_filt_group = merged_melt_filt_group_V08
#var_to_diff='kmeans_grouping'



library(dplyr)

wilcox_all_vars<-function(merged_melt_filt_group, var_to_diff){
  
  wilcox_stats_group<-merged_melt_filt_group %>%
    group_by(symbol) %>%
    do(w=wilcox.test(value~(!!sym(var_to_diff)), data=.) )%>%
    dplyr::summarize(symbol, Wilcox=w$p.value) %>% 
    mutate(p.adj=p.adjust(Wilcox)) %>%
    dplyr::filter(Wilcox<0.05) %>%
    dplyr::filter(p.adj<0.05) %>%
    arrange(Wilcox, decreasing=FALSE) %>%
    as.data.frame() 
  
  
  sum_meds<-merged_melt_filt_group %>%
   # group_by(symbol, kmeans_grouping) %>%
    dplyr::filter(symbol %in% wilcox_stats_group$symbol)%>%
   group_by(symbol, kmeans_grouping) %>%
    dplyr::summarise(median=median(value)) 
    
  

  
  sum_meds_dirs<-reshape2::dcast(sum_meds  ,symbol ~ kmeans_grouping, value.var = 'median') 
  sum_meds_dirs$direction=as.logical(sum_meds_dirs[2]> sum_meds_dirs[3])
   
  
  
  de_all_vars<-wilcox_stats_group[order(wilcox_stats_group$Wilcox),]
  de_all_vars$direction<-sum_meds_dirs[match( de_all_vars$symbol,sum_meds_dirs$symbol ),'direction' ]
  
  return(de_all_vars)
}


get_de_features_by_group<-function(merged_melt_filt_group, var_to_diff){
  #'
  #'
  #merged_melt_filt_group=merged_melt_filt_g1_ct
  var='kmeans_grouping'
  merged_melt_filt_group_V08<-merged_melt_filt_group%>%dplyr::filter(VISIT=='V08')
  g1_vs_control=wilcox_all_vars(merged_melt_filt_group_V08, var_to_diff = var_to_diff)
  return(g1_vs_control)
  
}


merged_melt_filt


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
    dplyr::summarize(symbol, Wilcox=w$p.value) %>% 
    mutate(p.adj=p.adjust(Wilcox)) %>%
    dplyr::filter(p.adj<0.05) %>%
    arrange(Wilcox, decreasing=FALSE) %>%
    as.data.frame() 
  
  
  most_sig_over_time_group<-wilcox_stats_group[order(wilcox_stats_group$Wilcox),]
  
  
  return(most_sig_over_time_group)
}




#se=se_rnas; feat_names=feats_rna_all  ; MOFAobject_clusts = MOFAobject_sel; 
#MOFAobject_sel
#se;feat_names
create_merged_melt<-function(se, feat_names, view, MOFAobject_clusts, run_cpm=TRUE){
  #' merged_melt_filt: is a merged sumamrized experiment, in long format, 
  #' that includes multiple visits
  #' and all molecules per patient 
  #' @param se
  #' @param
  #'  
  merged_melt_orig_1<-create_visits_df(se, clinvars_to_add, feat_names = feat_names, filter_common = TRUE, run_cpm=run_cpm)
  merged_melt_orig<-merged_melt_orig_1
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
    #' @merged
    merged_melt=merged_melt_orig[merged_melt_orig$COHORT==sel_cohort, ]
  }else{
    merged_melt=merged_melt_orig
  }
  merged_melt_pd=merged_melt_orig[merged_melt_orig$COHORT==1, ]
  merged_melt_ct=merged_melt_orig[merged_melt_orig$COHORT==2, ]
  
  merged_melt_pd<-merged_melt
  y_clust='NP2PTOT_LOG'
  clust_name= paste0(y_clust, '_clust')
  
  
 
  merged_melt$kmeans_grouping<-groups_from_mofa_factors(patnos=merged_melt$PATNO, MOFAobject_clusts= MOFAobject_clusts, y_clust)
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
  
  return(merged_melt_filt)
  
}
