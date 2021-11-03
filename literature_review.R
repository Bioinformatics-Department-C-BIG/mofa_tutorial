# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics

colname<-'Data'


level1<-c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Lipidomics', 'Metagenomics', 'miRNAs')


preprocessing<-function(stats_filter,colname){
  #' Split the column 
  #'
  
  omics_data<-str_split(stats_filter[[colname]], ',|\r|\n')
  omics_data<-lapply(omics_data,trimws)
  omics_data<-omics_data[!is.na(omics_data)]
  omics_data<-unlist(omics_data)
  omics_data<-omics_data[omics_data!='']
  
  
}



preprocessing_combinations<-function(stats_filter,x){
  #' Split the column 
  #'
  
  omics_data<-str_split(stats_filter[[colname]], ',|\r|\n\ ')
  omics_data<-lapply(omics_data,trimws)
  omics_data<-omics_data[!is.na(omics_data)]
  omics_data<-lapply(omics_data,unlist)
  omics_data<-omics_data[omics_data!='']
  
  #omics_data<-omics_data[!is.na(omics_data)]
  
  
  combinations<-lapply(omics_data, function(x) {
    x<-unlist(x)
    x<-x[x %in% level1]
    if (length(x)>2){
      x<-x[order(x)]
    combn(x,2, FUN=paste, collapse=' - ')}
    
  }
  
  )
  return(combinations)
  }


get_frequencies<-function(omics_data){
  
  #' Get the frequency of each occurrence 
  omics_data_frequencies<-table(tolower(unlist(omics_data)))
  omics_data_frequencies<-omics_data_frequencies[order(-omics_data_frequencies)]
  omics_data_frequencies<-data.frame(omics_data_frequencies)
  return(omics_data_frequencies)
  
  
}





library('readxl')
library('stringr')
library(ggplot2)
library(data.table)



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )

#stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_cancer_literature_curated.xlsx')

#write.csv(table(tolower(stats$Disease)),'Frequency_stats.csv')
stats_summarize<-as.data.frame(table(tolower(stats$Disease)))
ordered_stats<-stats_summarize[order(-stats_summarize$Freq),]
stats_summarize[order(stats_summarize$Var1),]


#write.csv(ordered_stats,'Frequency_stats_cancer.csv')

write.csv(ordered_stats,'Frequency_stats_not_cancer.csv')

ordered_stats
stats$`Objective-Code`

###Filters
#### 1. remove same sample 
#stats_filter<-stats[stats$`Data` %like% 'Proteomics',]


#### Split the omics and count number each used 
stats_filter<-stats[stats$`Same sample`=='Yes',]
omics_data<-preprocessing(stats_filter, 'Data')

omics_data_frequencies_1<-get_frequencies(omics_data)
omics_data_frequencies_1$Same_sample<-'Yes'

stats_filter<-stats[stats$`Same sample`=='No',]

omics_data<-preprocessing(stats_filter, 'Data')
omics_data_frequencies_2<-get_frequencies(omics_data)
#omics_data_frequencies_2<-omics_data_frequencies[omics_data_frequencies$Freq>1,]
omics_data_frequencies_2$Same_sample<-'No'

omics_data_frequencies<-rbind(omics_data_frequencies_2,omics_data_frequencies_1)
omics_data_frequencies<-omics_data_frequencies[order(-omics_data_frequencies$Freq),]

main_data2 <- omics_data_frequencies[ omics_data_frequencies$Var1 %in% tolower(level1), ]

omics_data_to_plot<-omics_data_frequencies


#x<-'Data'
#str_split(stats_filter[[x]], ',|\r|\n|\ ')
#barplot(omics_data_frequencies$Freq)
#omics_data_frequencies<-omics_data_frequencies[-1,]
ggplot(omics_data_to_plot, aes(x=Var1, y=Freq, fill=Same_sample))+
  geom_bar(stat='identity',position='stack')+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 30, vjust = 0.5, hjust=1))
  


####### Objectives 

stats_filter<-stats[stats$`Same sample`=='Yes',]


omics_data<-preprocessing(stats_filter, 'Objective-Code')


#### Group objectives more 
omics_data<-sapply(omics_data,function(x) gsub('.*Diagnosis.*|.*diagnosis.*', 'Diagnosis', x))


omics_data_frequencies_1<-get_frequencies(omics_data)
omics_data_frequencies_1$Same_sample<-'Yes'



#### Get combinations 
###Co-Occurence#

stats_filter<-stats[stats$`Same sample`=='Yes',]
omics_data_combinations<-preprocessing_combinations(stats_filter, 'Data')
omics_data_frequencies_1<-get_frequencies(omics_data_combinations)
omics_data_frequencies_1$Same_sample<-'Yes'
omics_data_frequencies_1<-omics_data_frequencies_1[omics_data_frequencies_1$Var1!=' - ',]



stats_filter<-stats[stats$`Same sample`=='No',]
omics_data_combinations<-preprocessing_combinations(stats_filter, 'Data')
omics_data_frequencies_2<-get_frequencies(omics_data_combinations)
omics_data_frequencies_2$Same_sample<-'No'

edge_list<-data.frame(do.call(rbind, str_split(omics_data_frequencies_1$Var1, ' - ')))
edge_list$weight<-omics_data_frequencies_1$Freq



    
ggplot(omics_data_frequencies_1, aes(x=Var1, y=Freq, fill=Same_sample))+
  geom_bar(stat='identity',position='stack')+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 20, vjust = 0.5, hjust=1))

library(igraph)
#edge_list<-edge_list[edge_list$weight>1,]
net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)



plot.igraph(net, edge.width=edge_list$weight)



