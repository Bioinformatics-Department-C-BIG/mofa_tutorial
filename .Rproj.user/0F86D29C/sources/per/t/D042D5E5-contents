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
  
  
  splitted<-str_split(df[[colname]], ',|\r|\n') # split by space, comma, newline
  splitted<-lapply(splitted,trimws)           # remove whitespace
  splitted<-splitted[!is.na(splitted)]      # remove nas 
  
  
  return(splitted)
  
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
stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )

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


### todo remove reviews? #
## GLOBAL FILTER
stats <- stats %>%
  filter(Type!= 'Review')



stats$same_sample<-as.factor(tolower(stats$same_sample))
  get_frequencies_by_group<-function(stats,colname){

  df_by_group<- stats %>%
    #filter(Type!= 'Review') %>%
    group_by(same_sample) %>%
    group_map(~ preprocessing(.x, colname) %>%
                get_frequencies() 
    )  %>%
    map_df(I, .id='same_sample')
  
    
    df_by_group$same_sample<-as.factor(df_by_group$same_sample)
    #### Filter by categories 
  
    df_by_group$same_sample<-recode_factor(df_by_group$same_sample,
                  '1'=cats[1], '2'=cats[2], '3'=cats[3], '4'=cats[4], '5'=cats[5] )

  return(df_by_group)
}

stats$same_sample<-as.factor(tolower(stats$same_sample))



colname<-'Data'
#### Also filter by omics
frequencies_by_group<-get_frequencies_by_group(stats, colname)

freq_to_plot<-frequencies_by_group %>% filter(Var1 %in% tolower(level1))


colname<-'Objective-Code'
frequencies_by_group<-get_frequencies_by_group(stats, colname)
frequencies_by_group$Freq<-as.numeric(frequencies_by_group$Freq)
frequencies_by_group$Var1<-sapply(frequencies_by_group$Var1,
                                  function(x) 
                                    gsub('.*diagnosis.*|*prognosis*', 'Diagnosis/Prognosis', x))

freq_to_plot<-frequencies_by_group



#### Create the stacked plot

library(grid)

ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL)+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,2,1.5),"cm"))

ggsave(paste0('plots/by', as.character(colname), '.png'), width = 8, height = 5)


#Plot by objective


# Plot by disease 


####### Objectives 

##### Plot by objective
#' 1. Group objectives to higher level


#### Get combinations 
### Co-Occurrence #

colname<-'Data'

# 
# get_combination_frequencies_by_group<-function(stats,colname){
#   df_by_group <- stats %>%
#     group_by(same_sample) %>%
#     group_map(~ preprocessing(.x, colname) %>%
#                 preprocessing_combinations()%>%
#                 get_frequencies() 
#     )  %>%
#    map_df(I, .id='same_sample')
#   
#   df_by_group$same_sample<-as.factor(df_by_group$same_sample)
#   #### Filter by categories 
#   
#   df_by_group$same_sample<-recode_factor(df_by_group$same_sample,
#                                          '1'=cats[1], '2'=cats[2], '3'=cats[3], '4'=cats[4], '5'=cats[5] )
#   
#   return(frequencies_by_group)
# }
# 
# 
# 
# comb_frequencies_by_group<-get_combination_frequencies_by_group(stats, 'Data')


omics_data<-df_by_group[[1]]
preprocessing_combinations<-function(omics_data){
  
  #' Create combinations of omics datasets  
  
  
  omics_data<-omics_data[omics_data!='']
  # 
  omics_data<-omics_data[!is.na(omics_data)]
  #' Create pairs of omics 
  combinations<-lapply(omics_data, function(x) {
    x<-unlist(x)
    x<-x[tolower(x) %in% tolower(level1)]
    if (length(x)>1){
      x<-x[order(x)]
      combn(x,2, FUN=paste, collapse=' - ')}
    
  }
  
  )
  return(unlist(combinations))
}



df_by_group <- stats %>%
  group_by(same_sample) %>%
  group_map(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations %>%
              get_frequencies() 
  )  %>%
  map_df(I, .id='same_sample')

df_by_group$same_sample<-as.factor(df_by_group$same_sample)
#### Filter by categories 

#df_by_group$same_sample<-recode_factor(df_by_group$same_sample,
 #                                      '1'=cats[1], '2'=cats[2], '3'=cats[3], '4'=cats[4], '5'=cats[5] )

comb_frequencies_by_group<-df_by_group

comb_frequencies_by_group<- comb_frequencies_by_group %>%
  group_by(same_sample) %>%
 filter(Freq>1)
    
comb_frequencies_by_group %>%
#filter(same_sample %in% c(2)) %>%
ggplot(data=., aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL, title='Combinations with > 1 occurences')+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 35, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,1.7,2.5),"cm"))

  


##### 
#' Create an edge list from the frequency table ! 
omics_data_frequencies<-comb_frequencies_by_group %>% filter(same_sample %in% c(2))
# aggregated or separately? 
aggr_freqs<-aggregate(Freq ~ Var1, omics_data_frequencies, sum)

comb_freq<-aggr_freqs

#comb_freq<- comb_frequencies_by_group %>% filter(same_sample ==3)
edge_list<-data.frame(do.call(rbind, str_split(comb_freq$Var1, ' - ')))
edge_list$weight<-comb_freq$Freq
edge_list<-edge_list[order(edge_list$weight, decreasing = TRUE),]


edge_list

library(igraph)
#edge_list<-edge_list[edge_list$weight>1,]
net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)

plot.igraph(net, edge.width=edge_list$weight)



####

#whay objectives  do proteomics and transcriptomics have? 
grepl('proteomics.*transcriptomics',tolower(stats$Data[23]))
stats$Data<-tolower(stats$Data)

new<-dplyr::filter(stats, grepl("proteomics.*transcriptomics",tolower(Data)))
new$`Objective-Code`

 
