library(dplyr)


metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)


combined$COHORT_DEFINITION

combined[which(combined$NHY==101),]$NHY<-NA



#### 
library(ggplot2)
combined_choose<-combined[combined$COHORT %in% c(1,2,4),]
combined_choose<-combined[combined$COHORT %in% c(1,2) & combined$INEXPAGE %in% c('INEXPD', 'INEXHC', 'INEXSNCA', 'INEXLRRK2'),]

combined$INEXPAGE
group_by_var<-'INEXPAGE'
ggplot(combined_choose)+
  geom_density(aes_string(x='AGE_AT_VISIT', fill=group_by_var,
                   color=group_by_var, group=group_by_var), 
               alpha=0.5)



ggplot(combined)+
  geom_histogram(aes(x=AGE_AT_VISIT, fill=as.factor(COHORT), color=COHORT))



library(ggplot2)

graphics.off()

### Create new variables from the averages 
## scopa does not have a total - maybe just add scopa total ? 
# because the average is the same as np3total

get_averages<-function(combined,sub_pattern ){
  # groups and averages specific columns 
    # TODO somehwte it is considering NAs as zeros CHECK 
    sub_pattern=paste0(sub_pattern,'[1-9]')  #For testing
    #sub_pattern='SCAU[1-9]'
  
    df<-combined[ , grepl( sub_pattern, colnames( combined ) )
                          & !grepl('TOT',  colnames( combined ) ) ]

    print(colnames(df))
    df<-as.data.frame(apply(df, 2, as.numeric))
    
    ind <- rowSums(is.na(df)) == ncol(df)

    df$sca_tot<-rowSums(df, na.rm = TRUE)
    # If all rows are NA then replace zero sum by NA
    #TODO: or do not run in the first place
    df$sca_tot[ind]<-NA 
    #combined$PATNO = as.factor(combined$PATNO )
    return(df$sca_tot)
    print(df$sca_tot)

}


sca_assess<-combined[ , grepl( sub_pattern, colnames( combined ) ) ]
colnames(sca_assess)
library(stringr)
sub_patterns=c( 'SCAU', 'STAIAD', 'NP3','NP1', 'NP2')
sub_pattern=paste0(sub_patterns[1],'[1-9]')  #For testing
sub_patterns_all<-paste(sub_patterns, collapse='|')
sub_patterns_all  
df$SCAU26C
df<-combined[ , grepl( sub_patterns_all, colnames( combined ) )
                & !grepl('TOT',  colnames( combined ) ) ]
  df=as.data.frame(df)
  df
  df=df[sapply(df, is.numeric)]
  df_log<-sapply(df, function(x) log2(x+10^-6))
  
  df_log
  df_log=as.data.frame(df_log)
  df_lognan_r<-sapply(df_log, is.nan)
  df_log[df_lognan_r]<-NA
  #combined_new<-merge(combined,df_log)
  combined_new<-mutate(df_log,combined)
  dim(combined_new)
  metadata_output_all<-paste0(output_files, 'combined_log',  '.csv')
  write.csv2(combined_new,metadata_output_all, row.names = FALSE)

  hist(combined_new$SCAU2)
  hist(combined$SCAU2)
  
  hist(df_log$SCAU2)
  
  # conevrt to apply 
#combined[,sub_pattern]<-
  
# ADD THE NEW averages   
avs<-sapply(sub_patterns,get_averages, combined=combined)
combined<-cbind(combined,avs)




#combined %>% 

  

scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU')

#

# postural instability gait disorder dominant
#### Create averages 
##




scales_in_stage<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU', 'STAIAD')




##### INSPECT DISTIRBUTIONS to decide on appropriate transforms 
###  Create plots only for the specific datasets
### If we are going to normalize/standardize by min-max it is better to do it only for the specific samples? 
i=2
graphics.off()
combined_p<-combined[combined$PATNO_EVENT_ID %in% common_samples[1:100], ]

#### Histograms to check the distributions of the clinical variables before and after processing 
for (i in 1:length(scales_in_stage)){
  
    ### 
  
    y=scales_in_stage[i]
    
  
    vals<-combined_p[,y]  
    hist(vals)
    
    vals_zero<-vals
    vals_zero[which(vals==0)]<-10^-6
    
    vals_test<-vals[!(is.na(vals) | vals==0)]
    #vals_test<-sample(vals_test, 5000)
    #hist(scale(vals_test, center = FALSE))
    #print(shapiro.test(vals_test))
    
    #hist(combined_to_plot[,y], na.rm=TRUE)
    
    png(paste0(outdir_orig,'metadata/hist_',y,'.png'))
    hist(vals, main=y)
    dev.off()
    
    png(paste0(outdir_orig,'metadata/hist_log_',y,'.png'))
    log_norm<-log2(vals)
    hist(log_norm, main=y)
    dev.off()
    
    log_norm<-log2(vals_zero)
    hist(log_norm)
    
    scaled<-scale(log2(vals_zero), center = FALSE)
    print(summary(scaled))
    hist(scaled)  
    
}


scales_in_stage<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU', 'STAIAD')



##### AVERAGE ALL THE DATA 
###
###

scaled=DataFrame();scaled_log<-DataFrame();scaled_log_sc<-DataFrame()

for (i in 1:length(scales_in_stage)){
  
  ### scale all the vars to average them 
  ### FOR NP3TP etc. make log 
  
  y=scales_in_stage[i]
  y2<-paste0(y,'_sc')
  comb_scale_data<-combined[,y]
  #print(paste(y, which(comb_scale_data==0)))
  
  comb_scale_data[which(comb_scale_data==0)]<-(10^-6)
  scaled[,y2]<-scale(comb_scale_data, center = FALSE)
  scaled_log[,y2]<-log2(comb_scale_data)
  
  scaled_log_sc[,y2]<-scale(log2(comb_scale_data), center=FALSE)
  
  
  

}

average_if_not_na<-function(df){
  ind <- rowSums(is.na(df)) == ncol(df)
  df$average<-rowMeans(as.data.frame(df), na.rm = TRUE)
  # If all rows are NA then replace zero sum by NA
  #TODO: or do not run in the first place
  df$average[ind]<-NA
  return(df$average)
}


#scaled$average=rowMeans(as.data.frame(scaled), na.rm=TRUE)

combined$STAGE_AV<-average_if_not_na(scaled)
combined$STAGE_LOG_AV<-average_if_not_na(as.data.frame(scaled_log))
combined$STAGE_LOG_SCALE_AV<-average_if_not_na(as.data.frame(scaled_log_sc))


hist(combined$STAGE_AV)



###### Now filter by relevant patients to make the plots 
### First filter by combined_filt\


combined_filt<-combined[combined$PATNO_EVENT_ID %in% common_samples[1:100], ]

combined_filt<-combined
# inspect patients
#View(combined_filt[combined_filt$PATNO=='3710',])


table(combined_filt$PAG_NAME_M3)
combined_filt$line_group = with(combined_filt,paste(PATNO,PAG_NAME_M3,PDSTATE, sep='_' ))


combined_filt$line_group
combined_filt$COHORT_DEFINITION

combined_filt$NHY
# FILTER OUT non visits
combined_filt<-combined_filt[grepl('V',combined_filt$EVENT_ID  ) | grepl('BL',combined_filt$EVENT_ID  ), ]

combined_filt$NHY
combined_filt=as.data.frame(combined_filt)


PS_101<-combined[which(combined$NHY==101),]


dim(PS_101[,c('COHORT_DEFINITION','NHY' )])

# REMOVE OUTLIERS FOR plot consistency
combined_filt<-combined_filt %>% 
            filter(NHY!=101) %>%
            filter(PAG_NAME_M3=='NUPDRS3')


#View(combined_filt[combined_filt$PATNO=='3710',])







outl<-max(combined_filt[y], na.rm=TRUE)
group='line_group'
x='EVENT_ID'
shape='PAG_NAME_M3'

#time points
scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU', 'STAGE_AV', 'STAGE_LOG_AV', 'STAGE_LOG_SCALE_AV')


tps<-read.csv(paste0('ppmi/ppmi_data/','visit_tps.csv'), sep=',')
tps#


tps[,1]<-gsub(' ','', tps[,1] )


##  
inds<-match(combined_filt$EVENT_ID, as.character( tps[,1]))
combined_filt$months<-as.numeric(tps[inds,2])
x='months'
combined$INEX
combined_to_plot<-combined_filt%>% select(c( y, x, group,
                                             scales, 'line_group', 'PATNO' , 'PAG_NAME_M3',
                                'AGE_AT_VISIT', 'COHORT_DEFINITION', 
                                'INEXPAGE', 'PDSTATE', 'PAG_NAME_M4'))
combined_filt$PD
combined_filt$COHORT_DEFINITION
combined_to_plot$COHORT_DEFINITION
fw<-'COHORT_DEFINITION'

combined_filt$line_group
### create an average of all clin vars 


combined_filt[combined_filt$COHORT==4,]$STAGE_LOG_SCALE_AV
group
y='STAGE_LOG_AV'
colour_by<-'PAG_NAME_M4'
colour_by<-'PAG_NAME_M3'
colour_by<-'PDSTATE'


scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU', 'STAGE_AV', 'STAGE_LOG_AV', 'STAGE_LOG_SCALE_AV')
scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCAU', 'STAGE_AV')

scales<-c( 'STAGE_AV', 'STAGE_LOG_AV', 'STAGE_LOG_SCALE_AV')
y='STAGE_AV'
formula_1<-as.formula('~COHORT_DEFINITION')
formula_1<-as.formula('~INEXPAGE')
formula_1<-as.formula('~PDSTATE')


combined_filt$PAG_NAME_M3
combined_to_plot<-combined_to_plot[combined_to_plot$INEXPAGE %in% c('INEXHC', 'INEXPD'),]

for (y in scales){

  p<-ggplot(combined_to_plot, aes_string( x=x, color=colour_by, group='line_group'))+
  geom_point(aes_string(y=y,color=colour_by, shape=shape))+
    geom_line(aes_string(y=y,color=colour_by, group=group))+
    guides( shape='none', group='none')#+
  
  
  
  p
  #theme(legend.position="none")
      #theme(legend.position="bottom", legend.text=element_text(size=2))+
    #theme(plot.margin=unit(c(-0.5, 1, 10, 0.5), units="line"))
          
  p+facet_wrap(formula_1, nrow = 4)
  p
  ggsave(paste0(outdir_orig,'metadata/lines_',paste0(formula_1, collapse=''),group, colour_by, y,'.jpeg' ), width=10, height=7)
  
  
  
  p<-ggplot(combined_to_plot, aes_string( x=x, color=colour_by, group='line_group'))+
    geom_point(aes_string(y=y,color=group, shape=shape))+
    #geom_line(aes_string(y=y,col=group, group=group)) +
    #geom_boxplot(aes_string(x=x, y=y))+
    geom_violin(aes_string(x=x, y=y, color=x))+
    
    theme(legend.position="none")
  p
  #theme(legend.position="bottom", legend.text=element_text(size=2))+
  #theme(plot.margin=unit(c(-0.5, 1, 10, 0.5), units="line"))
  
  p+facet_wrap(formula_1, nrow = 4)
  ggsave(paste0(outdir_orig,'metadata/box_',paste0(formula_1, collapse=''), y,'.jpeg' ), width=10, height=7)
  
  
  
  }

graphics.off()
  #geom_smooth(aes_string())
 # scale_y_continuous(limits = c(0, 7))



#### Create averages of all stages! 
# 

#### WHICH PRODROMAL ARE IN STAGE 3 
table(combined_filt[combined_filt$COHORT==4, 'STAGE_AV'])

combined_filt[combined_filt$COHORT==4 & combined_filt$STAGE_AV >1,]$PATNO


