

VISIT='V08'
metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)


combined$COHORT_DEFINITION

library(ggplot2)
### Create new variables from the averages 

get_averages<-function(combined,sub_pattern ){
  # groups and averages specific coluymns 
    sca_assess<-combined[ , grepl( sub_pattern, colnames( combined ) )
                          & !grepl('TOT',  colnames( combined ) ) ]

    sca_assess<-as.data.frame(apply(sca_assess, 2, as.numeric))
    sca_assess$sca_average<-rowMeans(sca_assess, na.rm = TRUE)
    #combined$PATNO = as.factor(combined$PATNO )

}
sub_pattern='NP1'  
sca_assess<-combined[ , grepl( sub_pattern, colnames( combined ) ) ]
colnames(sca_assess)

sub_patterns=c('NP1','NP3', 'NP2', 'NP4', 'SCA')


for (sub_pattern in sub_patterns){
  combined[,sub_pattern]<-get_averages(combined, sub_pattern)
  
}

#combined %>% 

  



scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY', 'SCA', 'NP1', 'NP2' , 'NP3', 'NP4')

#
#View(combined[combined$PATNO=='3710',])


### First filter by combined_filt\


combined_filt<-combined[combined$PATNO %in% common_samples[1:100], ]


# inspect patients
#View(combined_filt[combined_filt$PATNO=='3710',])


table(combined_filt$PAG_NAME_M3)
combined_filt$line_group = with(combined_filt,paste(PATNO,PAG_NAME_M3,PDSTATE, sep='_' ))
combined_filt$line_group
combined_filt$COHORT_DEFINITION


# FILTER OUT non visits
combined_filt<-combined_filt[grepl('V',combined_filt$EVENT_ID  ) | grepl('BL',combined_filt$EVENT_ID  ), ]


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
for (y in scales){

  p<-ggplot(combined_to_plot, aes_string( x=x, color='line_group', group='line_group'))+
  geom_point(aes_string(y=y,color=group, shape=shape))+
    geom_line(aes_string(y=y,col=group, group=group))+
      theme(legend.position="none")
  p+facet_wrap(~COHORT_DEFINITION, nrow = 3)
  ggsave(paste0(outdir_orig,'metadata/',y,'.jpeg' ), width=10, height=7)
}
  #geom_smooth(aes_string())
 # scale_y_continuous(limits = c(0, 7))
         
  

