library(dplyr)

VISIT='V08'

com
metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)


combined$COHORT_DEFINITION

combined[which(combined$NHY==101),]$NHY<-NA


library(ggplot2)
### Create new variables from the averages 

get_averages<-function(combined,sub_pattern ){
  # groups and averages specific coluymns 
    sca_assess<-combined[ , grepl( sub_pattern, colnames( combined ) )
                          & !grepl('TOT',  colnames( combined ) ) ]

    print(colnames(sca_assess))
    sca_assess<-as.data.frame(apply(sca_assess, 2, as.numeric))
    sca_assess$sca_average<-rowMeans(sca_assess, na.rm = TRUE)
    #combined$PATNO = as.factor(combined$PATNO )

}
sub_pattern='NP1'  
sca_assess<-combined[ , grepl( sub_pattern, colnames( combined ) ) ]
colnames(sca_assess)

sub_patterns=c('NP1','NP3', 'NP2', 'NP4', 'SCA', 'STAIAD')

#add gait 
for (sub_pattern in sub_patterns){
  combined[,sub_pattern]<-get_averages(combined, sub_pattern)
  
}

combined$STAIAD
#combined %>% 

  

scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCA', 'NP1', 'NP2' , 'NP3', 'NP4')

#

# postural instability gait disorder dominant
#### Create averages 
##




scales_in_stage<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCA', 'STAIAD')



i=2
graphics.off()
for (i in 1:length(scales_in_stage)){
  
    ### 
  
    y=scales_in_stage[i]
    
  
    vals<-combined[,y]  
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


scales_in_stage<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCA', 'STAIAD')


### AVERAGE ALL THE DATA 
scaled=DataFrame();scaled_log<-DataFrame();scaled_log_sc<-DataFrame()

for (i in 1:length(scales_in_stage)){
  
  ### scale all the vars to average them 
  ### FOR NP3TP etc. make log 
  
  y=scales_in_stage[i]
  y2<-paste0(y,'_sc')
  comb_scale_data<-combined[,y]
  #print(paste(y, which(comb_scale_data==0)))
  
  comb_scale_data[which(comb_scale_data==0)]<(10^-6)
  scaled[,y2]<-scale(comb_scale_data, center = FALSE)
  scaled_log[,y2]<-log2(comb_scale_data)
  
  scaled_log_sc[,y2]<-scale(log2(comb_scale_data), center=FALSE)
  
  
  

}

scaled$average=rowMeans(as.data.frame(scaled), na.rm=TRUE)
combined$STAGE_AV<-scaled$average
combined$STAGE_LOG_AV<-rowMeans(as.data.frame(scaled_log), na.rm=TRUE)
combined$STAGE_LOG_SCALE_AV<-rowMeans(as.data.frame(scaled_log_sc), na.rm=TRUE)


hist(combined$STAGE_AV)


#View(combined[combined$PATNO=='3710',])


### First filter by combined_filt\


combined_filt<-combined[combined$PATNO %in% common_samples[1:100], ]


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


dim(PS_101[,c('COHORT_DEFINITION','NHY' , '')])

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
scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCA', 'NP1', 'NP2' , 'NP3', 'NP4', 'STAGE_AV', 'STAGE_LOG_AV', 'STAGE_LOG_SCALE_AV')


tps<-read.csv(paste0('ppmi/ppmi_data/','visit_tps.csv'), sep=',')
tps#


tps[,1]<-gsub(' ','', tps[,1] )


##  
inds<-match(combined_filt$EVENT_ID, as.character( tps[,1]))
combined_filt$months<-as.numeric(tps[inds,2])
x='months'

combined_to_plot<-combined_filt%>% select(c( y, x, group,
                                             scales, 'line_group', 'PATNO' , 'PAG_NAME_M3',
                                'AGE_AT_VISIT', 'COHORT_DEFINITION'))


combined_filt$COHORT_DEFINITION
combined_to_plot$COHORT_DEFINITION
fw<-'COHORT_DEFINITION'


### create an average of all clin vars 


combined_filt[combined_filt$COHORT==4,]$STAGE_LOG_SCALE_AV

y='STAGE_LOG_AV'

for (y in scales){

  p<-ggplot(combined_to_plot, aes_string( x=x, color='line_group', group='line_group'))+
  geom_point(aes_string(y=y,color=group, shape=shape))+
    geom_line(aes_string(y=y,col=group, group=group)) +
  theme(legend.position="none")
      #theme(legend.position="bottom", legend.text=element_text(size=2))+
    #theme(plot.margin=unit(c(-0.5, 1, 10, 0.5), units="line"))
          
  p+facet_wrap(~COHORT_DEFINITION, nrow = 3)
  ggsave(paste0(outdir_orig,'metadata/lines_',y,'.jpeg' ), width=10, height=7)

  
  
  p<-ggplot(combined_to_plot, aes_string( x=x, color='line_group', group='line_group'))+
    geom_point(aes_string(y=y,color=group, shape=shape))+
    #geom_line(aes_string(y=y,col=group, group=group)) +
    #geom_boxplot(aes_string(x=x, y=y))+
    geom_violin(aes_string(x=x, y=y, color=x))+
    
    theme(legend.position="none")
  p
  #theme(legend.position="bottom", legend.text=element_text(size=2))+
  #theme(plot.margin=unit(c(-0.5, 1, 10, 0.5), units="line"))
  
  p+facet_wrap(~COHORT_DEFINITION, nrow = 3)
  ggsave(paste0(outdir_orig,'metadata/box_',y,'.jpeg' ), width=10, height=7)
  
  
  
  }

graphics.off()
  #geom_smooth(aes_string())
 # scale_y_continuous(limits = c(0, 7))



#### Create averages of all stages! 
# 

#### WHICH PRODROMAL ARE IN STAGE 3 
table(combined_filt[combined_filt$COHORT==4, 'STAGE_AV'])

combined_filt[combined_filt$COHORT==4 & combined_filt$STAGE_AV >1,]$PATNO


