# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics



install.packages("readxl")
library('readxl')
stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer.xlsx' )

#stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_cancer_literature_curated.xlsx')

#write.csv(table(tolower(stats$Disease)),'Frequency_stats.csv')
stats_summarize<-as.data.frame(table(tolower(stats$Disease)))
ordered_stats<-stats_summarize[order(-stats_summarize$Freq),]
stats_summarize[order(stats_summarize$Var1),]


#write.csv(ordered_stats,'Frequency_stats_cancer.csv')

write.csv(ordered_stats,'Frequency_stats_not_cancer.csv')
