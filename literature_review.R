# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics

colname<-'Data'
library('dplyr')
library('purrr')
level1<-c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Lipidomics', 'Metagenomics', 'miRNAs')

# Process; if methylation or histone; add epigenomics!
preprocessing<-function(df,colname){
  #' Split the column 
  #' Return the split variables to get the frequencies 
  
  
  omics_data<-str_split(df[[colname]], ',|\r|\n') # split by space, comma, newline
  omics_data<-lapply(omics_data,trimws)           # remove whitespace
  omics_data<-omics_data[!is.na(omics_data)]      # remove nas 
  
  
  return(omics_data)
  
}

x
#omics_data<-frequencies_by_group[[1]]
x<-omics_data[[1]]
preprocessing_combinations<-function(omics_data){
  #' Split the column 
  #'
  
  #
  # omics_data<-str_split(stats_filter[[colname]], ',|\r|\n\ ')
  # omics_data<-lapply(omics_data,trimws)
  # omics_data<-omics_data[!is.na(omics_data)]
  
  #omics_data<-lapply(omics_data,unlist)
  omics_data<-omics_data[omics_data!='']
  # 
  omics_data<-omics_data[!is.na(omics_data)]
  #' Create pairs of omics 
  combinations<-lapply(omics_data, function(x) {
    x<-unlist(x)
    x<-x[tolower(x) %in% tolower(level1)]
    if (length(x)>2){
      x<-x[order(x)]
    combn(x,2, FUN=paste, collapse=' - ')}
    
  }
  
  )
  print(combinations)
  return(combinations)
  }

comb_frequencies_by_group<-get_combination_frequencies_by_group(stats, 'Data')


get_frequencies<-function(omics_data){
  
  #' Get the frequency of each occurrence 
  #' 
  #' 
  
  omics_data<-unlist(omics_data) # collapse accross studies 
  omics_data<-omics_data[omics_data!='']
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
#' 1. Split to same sample yes or no
#' 2 Simplify to lower level 
#' 3. Plot 
#' 
#' 
#' 
#' 
#' 


get_frequencies_by_group<-function(stats,colname){
  frequencies_by_group<- stats %>%
    group_by(`Same sample`) %>%
    group_map(~ preprocessing(.x, colname) %>%
                get_frequencies() 
    )  %>%
    map_df(I, .id='Same_sample')
  
  
    selected_cats<-c('yes', 'no', NA)
    categories<-levels(stats$`Same sample`)
    
    
    #### Filter by categories 
    sel<-which(categories %in% selected_cats)
    frequencies_by_group<-frequencies_by_group %>% filter(Same_sample %in% sel)
    
  return(frequencies_by_group)
}
stats$`Same sample`<-as.factor(tolower(stats$`Same sample`))



colname<-'Data'


frequencies_by_group<-get_frequencies_by_group(stats, colname)


#### Also filter by omics
freq_to_plot<-frequencies_by_group %>% filter(Var1 %in% tolower(level1))


#### Create the stacked plot

library(grid)

ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=Same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL)+
  scale_fill_discrete(labels=categories[sel])+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,1.7,1.2),"cm"))

ggsave('plots/byData.png', width = 8, height = 5)


#Plot by objective


# Plot by disease 
colname<-'Objective-Code'
frequencies_by_group<- stats %>%
  group_by(`Same sample`) %>%
  group_map(~ preprocessing(.x, colname) %>%
              get_frequencies() 
  )  %>%
  map_df(I, .id='Same_sample')
frequencies_by_group<-as.data.frame(frequencies_by_group)
frequencies_by_group$Freq<-as.numeric(frequencies_by_group$Freq)
frequencies_by_group$Var1<-sapply(frequencies_by_group$Var1,function(x) gsub('.*diagnosis.*|*prognosis*', 'Diagnosis/Prognosis', x))


freq_filt<-frequencies_by_group %>% filter(Same_sample ==3)
ggplot(freq_filt, aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=Same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL)+
  scale_fill_discrete(labels=categories[sel])+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,1.7,1.2),"cm"))


####### Objectives 

##### Plot by objective
#' 1. Group objectives to higher level


stats_filter<-stats[stats$`Same sample`=='Yes',]


omics_data<-preprocessing(stats_filter, 'Objective-Code')

#### Group objectives more 


omics_data_frequencies_1<-get_frequencies(omics_data)
omics_data_frequencies_1$Same_sample<-'Yes'


ggplot(omics_data_frequencies_1, aes(x=Var1, y=Freq, fill=Same_sample))+
  geom_bar(stat='identity',position='stack')+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 30, vjust = 0.5, hjust=1))

ggsave('plots/Objectives.png')

#### Get combinations 
### Co-Occurrence #

colname<-'Data'


get_combination_frequencies_by_group<-function(stats,colname){
  frequencies_by_group<- stats %>%
    group_by(`Same sample`) %>%
    group_map(~ preprocessing(.x, colname) %>%
                preprocessing_combinations()%>%
                get_frequencies() 
    )  %>%
   map_df(I, .id='Same_sample')
  
  #selected_cats<-c('yes')
  #categories<-levels(stats$`Same sample`)
  
  
  #### Filter by categories 
  #sel<-which(categories %in% selected_cats)
  #frequencies_by_group<-frequencies_by_group %>% filter(Same_sample %in% sel)
  
  return(frequencies_by_group)
}



comb_frequencies_by_group<-get_combination_frequencies_by_group(stats, 'Data')

comb_frequencies_by_group

  comb_frequencies_by_group %>%
  group_by(Same_sample) %>%
 filter(Freq>1) %>%
 filter(Same_sample %in% c('3')) %>%
ggplot(data=., aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=Same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL, title='Combinations with > 1 occurences')+
  scale_fill_discrete(labels=categories[sel])+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 35, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,1.7,2.5),"cm"))

  

#DO NOT SPLIT#
stats_filter=stats
omics_data_combinations<-preprocessing_combinations(stats_filter, 'Data')
omics_data_frequencies<-get_frequencies(omics_data_combinations)
omics_data_frequencies<-omics_data_frequencies[omics_data_frequencies$Var1!=' - ',]


##### 
#' Create an edge list from the frequency table ! 
omics_data_frequencies<-comb_frequencies_by_group
# aggregated or separately? 
aggr_freqs<-aggregate(Freq ~ Var1, comb_frequencies_by_group, sum)

comb_freq<-aggr_freqs

#comb_freq<- comb_frequencies_by_group %>% filter(Same_sample ==3)
edge_list<-data.frame(do.call(rbind, str_split(comb_freq$Var1, ' - ')))
edge_list
edge_list$weight<-comb_freq$Freq
edge_list<-edge_list[order(edge_list$weight, decreasing = TRUE),]


edge_list

library(igraph)
#edge_list<-edge_list[edge_list$weight>1,]
net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)

plot.igraph(net, edge.width=edge_list$weight)






####

#whay objectivesdo proteomics and transcriotmics have? 
grepl('proteomics.*transcriptomics',tolower(stats$Data[23]))
stats$Data<-tolower(stats$Data)

new<-dplyr::filter(stats, grepl("proteomics.*transcriptomics",tolower(Data)))
new$`Objective-Code`

 
