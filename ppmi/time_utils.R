











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
  
  clip_outliers<-function(X, x_times=1.5, upper=TRUE, lower=TRUE){
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




