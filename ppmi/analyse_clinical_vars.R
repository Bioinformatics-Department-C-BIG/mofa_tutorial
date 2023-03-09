

VISIT='V08'
metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)

metadata_output_all<-paste0(output_files, 'combined',  '.csv')
combined<-read.csv2(metadata_output_all)


library(ggplot2)



#combined$PATNO = as.factor(combined$PATNO )




scales<-c('NP1RTOT','NP2PTOT' , 'NP3TOT', 'NP4TOT', 'NHY')
table(combined$PAG_NAME_M3)
View(combined[combined$PATNO=='3710',])


### First filter by combined_filt\


combined_filt<-combined[combined$PATNO %in% common_samples[1:100], ]


# inspect patients
View(combined_filt[combined_filt$PATNO=='3710',])


table(combined_filt$PAG_NAME_M3)
combined_filt$line_group = with(combined_filt,paste(PATNO,PAG_NAME_M3,PDSTATE, sep='_' ))
combined_filt$line_group



# FILTER OUT non visits
combined_filt<-combined_filt[grepl('V',combined_filt$EVENT_ID  ) | grepl('BL',combined_filt$EVENT_ID  ), ]


# REMOVE OUTLIERS FOR plot consistency
combined_filt<-combined_filt %>% 
  filter(NHY!=101) %>%
  filter(PAG_NAME_M3=='NUPDRS3')


#View(combined_filt[combined_filt$PATNO=='3710',])







y=scales[1]
outl<-max(combined_filt[y], na.rm=TRUE)
group='line_group'
x='EVENT_ID'
shape='PAG_NAME_M3'
combined_to_plot<-combined_filt%>% select(c( y, x, group,
                                          scales, 'line_group', 'PATNO' , 'PAG_NAME_M3'))







ggplot(combined_to_plot, aes_string( x=x, color='line_group', group='line_group'))+
geom_point(aes_string(y=y,color=group, shape=shape))+
  geom_line(aes_string(y=y,col=group, group=group))+
    theme(legend.position="none")

ggsave(paste0(outdir_orig,'metadata/',y,'.jpeg' ))
  #geom_smooth(aes_string())
 # scale_y_continuous(limits = c(0, 7))
         
  

