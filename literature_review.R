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



get_frequencies<-function(x){
  
  #' Get the frequency of each occurrence 
  #' 
  #' 
  
  x<-unlist(x) # collapse accross studies 
  x<-x[x!='']
  x<-table(tolower(unlist(x)))
  x<-data.frame(x)
  #return(omics_data_frequencies)
  
  
}





library('readxl')
library('stringr')
library(ggplot2)
library(data.table)



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )

#stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_cancer_literature_curated.xlsx')


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


group_objectives<-function(df, Var1){
  #'Group objective code column 
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     gsub('.*diagnosis.*|*prognosis*', 'Diagnosis/Prognosis', x))
  return(df)
}

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
  return(combinations)
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





####
#' Combinations and their objectives
#' 
#' 


get_combs<-function(x) {
  x=unlist(lapply(x,trimws) )
  x<-x[tolower(x) %in% tolower(level1)]
  #' Get the combinations of omnics in pairs of 2
  if (length(x)>1){
    x<-x[order(x)]
    res<-combn(x,2, FUN=paste, collapse=' - ')
    return(res)
  }
  
}

library(tidyverse)




########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 

new<-stats %>% 
  mutate(`Objective-Code`=strsplit(`Objective-Code`, ',|\r|\n' ))%>%
  unnest(`Objective-Code`) 


new['Objective-Code'] <-apply(new['Objective-Code'], 1, function(x) trimws(tolower(x)))
colnames(new)[which(colnames(new)=='Objective-Code')]<-'objective'
colname='Data'
new<-group_objectives(new, 'objective')

new$objective<-as.factor(new$objective)
cats<-levels(new$objective)


keys<-pull(new %>%
             group_by(objective) %>%
             group_keys())


#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 


df_by_group<-new %>%
  group_by(objective) %>%
  group_map(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations() %>%
              get_frequencies() 
  )  %>%
  map_df(I, .id='objective') 




df_by_group<-as.data.frame(as.matrix(df_by_group))

df_by_group$objective<-as.numeric(df_by_group$objective)

key_names<-c(keys[df_by_group$objective])

df_by_group<-cbind(key_names,df_by_group)
df_by_group$Var1<-as.factor(df_by_group$Var1)


df_by_group$Freq<-as.numeric(df_by_group$Freq)

df_to_plot<-df_by_group %>% 
  group_by(Var1)  %>% 
  filter( sum(Freq) >= 7) %>% 
  group_by(objective)  %>% 
  filter( sum(Freq) >= 4) 


df_to_plot<-df_to_plot[!is.na(df_to_plot$key_names),]

plotbyObjective(df_to_plot )




plotbyObjective<-function(df){ 
  ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=Freq, fill=Var1))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,2,1.7),"cm"))
  save(paste0('plots/by', as.character(colname), '.png'), width = 8, height = 5)
  
  
}

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

#what objectives  do proteomics and transcriptomics have? 
grepl('proteomics.*transcriptomics',tolower(stats$Data[23]))

stats$Data<-tolower(stats$Data)
pattern1<-'*proteomics*transcriptomics*|*transcriptomics.*proteomics*'
pattern1<-'*proteomics*transcriptomics*|*transcriptomics.*proteomics*'
pattern1<-'*metabolomics.*proteomics*|*proteomics*metabolomics.*'

new<-dplyr::filter(stats, grepl(pattern1,tolower(Data)))


colname<-'Objective-Code'
freq_to_plot<-new %>% 
  preprocessing(colname) %>%
  get_frequencies() %>% 
  group_objectives('Var1')



plotby<-function(freq_to_plot, colname){ 
  ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,2,1.5),"cm"))
  
}

plotby(freq_to_plot)
## Plot the graph of proteomics transc
freq_to_plot







