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



group_objectives<-function(df, Var1){
  #'Group objective code column 
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     gsub('.*diagnosis.*|*prognosis*', 'Diagnosis/Prognosis', x))
  return(df)
}



library('readxl')
library('stringr')
library(ggplot2)
library(data.table)



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
stats<-read_excel('H:/My Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )


###Filters
#### 1. remove same sample 
#stats_filter<-stats[stats$`Data` %like% 'Proteomics',]


#### Split the omics and count number each used 
#' 1. Split to same sample yes or no
#' 2 Simplify to lower level 
#' 3. Plot 
#' 
#' 


### Remove reviews, remove rejected articles
## GLOBAL FILTER
stats <- stats %>%
  filter(Type!= 'Review')%>%
  filter(is.na(`Rejection /Critic`))



stats$same_sample<-as.factor(tolower(stats$same_sample))
get_frequencies_by_group<-function(stats,colname){
  
  df_by_group<- stats %>%
    group_by(same_sample) %>%
    group_map(~ preprocessing(.x, colname) %>%
                get_frequencies() 
    )  %>%
    map_df(I, .id='same_sample')
  
  
  return(df_by_group)
}

stats$same_sample<-as.factor(tolower(stats$same_sample))



colname<-'Data'
#### Also filter by omics
frequencies_by_group<-get_frequencies_by_group(stats, colname)

freq_to_plot<-frequencies_by_group %>% filter(Var1 %in% tolower(level1))


#### Create the stacked plot

library(grid)


plotByData(freq_to_plot)

ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL)+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,2,1.5),"cm"))

ggsave(paste0('plots/SingleOmicsby', as.character(colname), '.png'), width = 8, height = 5)



####### Objectives 

##### Plot by objective
#' 1. Group objectives to higher level


#### Get combinations 
### Co-Occurrence #

colname<-'Data'

omics_data<-df_by_group[[1]]


####
#' Combinations and their objectives
#' 
#' 

get_combs<- function(x){
    x<-unlist(x)
    x<-x[tolower(x) %in% tolower(level1)]
    if (length(x)>1){
      x<-x[order(x)]
      combn(x,2, FUN=paste, collapse=' - ')}
    
  }

preprocessing_combinations<-function(x){
  
  #' Create combinations of omics datasets  
  
  
  x<-x[x!='']
  # 
  x<-x[!is.na(x)]
  #' Create pairs of omics 
  #' #
  #'
  combinations<-lapply(x,get_combs)
  return(combinations)
}



df_by_group <- stats %>%
  group_by(same_sample) %>%
  group_map(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations %>%
              get_frequencies() 
  )  %>%
  map_df(I, .id='same_sample')

freq_cutoff<-4

df_by_group<-df_by_group %>% 
  group_by(Var1)  %>% 
  filter( sum(Freq) >= freq_cutoff) 

plotByData(df_by_group)

plotByData<-function(df_by_group){
  ggplot(df_by_group, aes(x=reorder(Var1, -Freq, sum), y=Freq, fill=same_sample))+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL, title=paste0('Combinations with > ',freq_cutoff, ' occurences'))+
  theme(axis.text.x = element_text(size=rel(1.3),angle = 35, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,1.7,2.5),"cm"))
  ggsave(paste0('plots/byCombinations', as.character(colname), '.png'), width = 8, height=6)
  
}



library(tidyverse)

########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 
colnames(stats)[which(colnames(stats)=='Objective-Code')]<-'objective'
colnames(stats)[which(colnames(stats)=='Integration method-Category')]<-'method'

new<-stats %>% 
  mutate(objective=strsplit(objective, ',|\r|\n' ))%>%
  unnest(objective) 

new<-stats %>% 
  mutate(method=strsplit(method, ',|\r|\n' ))%>%
  unnest(method) 

x_group<-'objective'
x_group<-'method'

new<-stats %>% 
  mutate(objective=strsplit(objective, ',|\r|\n' ))%>%
  unnest(objective) 

colname='Data'

new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))


new<-group_objectives(new, x_group)

keys<-pull(new %>%
             group_by_at(x_group) %>%
             group_keys())


#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
df_by_group<-new %>%
  group_by_at(x_group) %>%
  group_map(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations() %>%
              get_frequencies() 
  )  %>%
  map_df(I, .id=x_group) 



# Attach the key names back to the dataframe 
df_by_group<-as.data.frame(as.matrix(df_by_group))
df_by_group[,x_group]<-as.numeric(df_by_group[,x_group])
key_names<-c(keys[df_by_group[,x_group]])
df_by_group<-cbind(key_names,df_by_group)


df_by_group$Freq<-as.numeric(df_by_group$Freq)

df_to_plot<-df_by_group %>% 
  group_by(Var1)  %>% 
  filter( sum(Freq) >= 8) %>% 
  group_by_at(x_group)  %>% 
  filter( sum(Freq) >= 4) 


df_to_plot<-df_to_plot[!is.na(df_to_plot$key_names),]

show_p<-plotbyObjective(df_to_plot )


show_p

plotbyObjective<-function(df){ 
  ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=Freq, fill=Var1))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,2,1.7),"cm"))
  
  ggsave(paste0('plots/by', as.character(colname), '.png'), width = 8, height=6)
  
  
  
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
edge_list<-edge_list[order(edge_l------------------------------------------+
                             ist$weight, decreasing = TRUE),]


edge_list

library(igraph)


net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)
save(paste0('plots/network', as.character(colname), '.png'), width = 8, height=6)

plot.igraph(net, edge.width=edge_list$weight)













