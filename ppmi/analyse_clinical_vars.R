

VISIT='V08'
metadata_output<-paste0(output_files, 'combined_', VISIT,  '.csv')
combined_bl<-read.csv2(metadata_output)


library(ggplot2)



combined$PATNO = as.factor(combined$PATNO )
ggplot(combined, aes(y=NP3TOT, x=EVENT_ID, fill=PATNO))+
  geom_point(aes(color=PATNO) )


combined$PATNO=as.factor(combined$PATNO)
ggplot(combined, aes(y=NP3TOT, x=EVENT_ID, fill=PATNO))+
  geom_line(aes(color=as.factor(PATNO) ))


