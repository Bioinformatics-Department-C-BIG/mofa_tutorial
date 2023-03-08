

VISIT='V08'
metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)


library(ggplot2)



#combined$PATNO = as.factor(combined$PATNO )

combined$PATNO=as.factor(combined$PATNO)
ggplot(combined, aes(y=NP3TOT, x=EVENT_ID, fill=PATNO))+
  geom_line(aes(color=as.factor(PATNO) ))



scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY')



### First filter by combined_filt\


combined_filt<-combined[combined$PATNO %in% common_samples[1:100], ]

combined_filt$REC_ID.x
combined_filt$line_group = with(combined_filt,paste(EVENT_ID,PATNO,PAG_NAME.x, PAG_NAME.y, sep='_' ))

y=scales[5]
outl<-max(combined_filt[y], na.rm=TRUE)
which[combined_filt[y]]
group='PATNO'
x='EVENT_ID'
inst=''
combined_filt<-combined_filt%>% select(c( y, x, group, scales, 'line_group' ))

combined_filt$PATNO=as.factor(combined_filt$PATNO)


combined_filt[c('NHY', x)]

combined_filt<-combined_filt[!(combined_filt$NHY==101),]


ggplot(combined_filt, aes_string( x=x, color='line_group', group='line_group'))+
geom_point(aes_string(y=y,color=group))+
  geom_line(aes_string(y=y,col=group, group=group))+
    theme(legend.position="none")
  #geom_smooth(aes_string())
 # scale_y_continuous(limits = c(0, 7))
         
  

