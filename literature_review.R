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


df<-new
Var1<-'objective'
group_objectives<-function(df, Var1){
  #'Group objective code column 
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     gsub('.*diagnosis.*|.*prognosis.*', 'Diagnosis/Prognosis', tolower(x)))
  return(df)
}


#TODO: make a function to check if there is methylomics
group_omics<-function(df, Var1){
  #'Group objective code column 
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     gsub('*methyl*', 'Epigenomics', tolower(x)))
  return(df)
}





library('readxl')
library('stringr')
library(ggplot2)
library(data.table)



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
stats<-read_excel('H:/My Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )
stats<-read_excel('E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )


###Filters
#### 1. remove same sample 
#stats_filter<-stats[stats$`Data` %like% 'Proteomics',]


#### Split the omics and count number each used 
#' 1. Split to same sample yes or no
#' 2 Simplify to lower level 
#' 3. Plot 
#' 
#' 

stats$Cancer<-c(rep('no',289), rep('yes',(nrow(stats)-289)))

### Remove reviews, remove rejected articles
## GLOBAL FILTER
stats <- stats %>%
  filter(Type!= 'Review')%>%
  filter(is.na(`Rejection /Critic`))%>%
  filter(tolower(same_sample)!='no')



stats$same_sample<-as.factor(tolower(stats$same_sample))
stats$Cancer<-as.factor(tolower(stats$Cancer))

x_group<-'Cancer'

get_frequencies_by_group<-function(stats,colname){
  
  df_by_group<- stats %>%
    group_by_at(x_group) %>%
    group_map(~ preprocessing(.x, colname) %>%
                get_frequencies() 
    )  %>%
    map_df(I, .id=x_group)
  
  
  return(df_by_group)
}

stats$same_sample<-as.factor(tolower(stats$same_sample))



colname<-'Data'
#### Also filter by omics
frequencies_by_group<-get_frequencies_by_group(stats, colname)

freq_to_plot<-frequencies_by_group %>% filter(Var1 %in% tolower(level1))
single_omics_frequencies=freq_to_plot

#### Create the stacked plot

library(grid)

plotByData<-function(df_by_group){
  ggplot(df_by_group, 
         aes(x=reorder(Var1, -Freq, sum), y=Freq))+
    aes_string(fill=y_group)+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL, title=paste0('Combinations with > ',freq_cutoff, ' occurences'))+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 35, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,1.7,2.5),"cm"))
  ggsave(paste0('plots/byCombinations', as.character(colname), '.png'), width = 8, height=6)
  
}




plotByData(freq_to_plot)

ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq))+
         aes_string(fill=x_group)+
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
    
  #' return omics combinations as individual strings to count them+
  x<-unlist(x)
    x<-x[tolower(x) %in% tolower(level1)]
    if (length(x)>1){
      x<-x[order(x)]
      combn(x,2, FUN=paste, collapse=' - ')
    }
}

preprocessing_combinations(preprocessing(stats, 'Data'))

preprocessing_combinations<-function(x){
  #' Create combinations of omics datasets  
  x<-x[x!='']
  x<-x[!is.na(x)]
  #' Create pairs of omics 
  #' #
  #'
  combinations<-lapply(x,get_combs)
  return(combinations)
}


total<-NROW(stats[!is.na(stats$Data),]$PMID)

y_group='same_sample'
y_group='Cancer'


df_by_group <- stats %>%
  #group_by(same_sample) %>%
  group_by_at(y_group) %>%
  group_map(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations %>%
              get_frequencies() 
  )  %>%
  map_df(I, .id=y_group)

freq_cutoff<-7

df_by_group_filtered<-df_by_group %>% 
  group_by(Var1)  %>% 
  filter( sum(Freq) >= freq_cutoff) 

plotByData(df_by_group_filtered)


combinations<-df_by_group


combinations <-df_by_group %>% separate(Var1, c("Omics1","Omics2"), sep = " - ")

ggplot(combinations)+aes(Omics1, Omics2, fill=abs(Freq)) +
  geom_tile()+
  geom_text(aes(label = round(Freq, 2)), size=7)+
  theme(axis.text.x = element_text(size=rel(1.5),angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text( size=rel(1.5)))
  

#TODO: add red find overflow 
  
  

comb_frequencies_by_group<-df_by_group


library(tidyverse)

########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 
colnames(stats)[which(colnames(stats)=='Objective-Code')]<-'objective'
colnames(stats)[which(colnames(stats)=='Integration method-Category')]<-'method'


#change here to select by objective or by method 
# TODO: MAKE this one variable to choose method or objective

x_group<-'method'
# and here!!! 
new<-stats %>% 
  mutate(method=strsplit(method, ',|\r|\n' ))%>%
  unnest(method) 


x_group<-'objective'
new<-stats %>% 
  mutate(objective=strsplit(objective, ',|\r|\n' ))%>%
  unnest(objective) 

colname='Data'

new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))
new<-group_objectives(new, 'objective')

keys<-pull(new %>%
             group_by_at(x_group) %>%
             group_keys())


#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
df_by_group<-new %>%
  group_by_at(x_group) 

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
filter( sum(Freq) >= 7) %>%
group_by_at(x_group)  %>%
filter( sum(Freq) >= 3)

#df_to_plot<-df_by_group
df_to_plot<-df_to_plot[!is.na(df_to_plot$key_names),]

plotbyObjective<-function(df){ 
  g<-ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=Freq, fill=Var1))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,2,1.7),"cm"))
  
  ggsave(paste0('plots/barplot_byGroup', as.character(x_group), '.png'), width = 8, height=6)
  return(g)
  
  
}


show_p<-plotbyObjective(df_to_plot )


show_p



##### 
#' Create an edge list from the frequency table ! 
# aggregated or separately? 


comb_freq<- comb_frequencies_by_group #%>% filter(Cancer ==2)
single_omics_frequencies_filtered<-single_omics_frequencies #%>% filter(Cancer ==2)

edge_list<-data.frame(do.call(rbind, str_split(comb_freq$Var1, ' - ')))
edge_list$weight<-comb_freq$Freq
edge_list<-edge_list[order(edge_list$weight, decreasing = TRUE),]


edge_list

library(igraph)

#### Add the frequency of single omics to the network vertices
df<-single_omics_frequencies_filtered

net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)
net_att<-df[match(V(net)$name, df$Var1),]


# TODO: assign omics frequencies of all samples or only cancer? 
vertex_attr(net, 'freq', index=V(net))<-single_omics_frequencies$Freq

p<-plot.igraph(net, edge.width=edge_list$weight, vertex.size=net_att$Freq)
save(p,paste0('plots/network', as.character(colname), '.png'))


g <- set.vertex.attribute(g,'id',1,'first_id')












