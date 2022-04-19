# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics
colname<-'Data'
library('dplyr')
library('purrr')
library(RColorBrewer)
source('utils.R')



### Load the package or install if not present
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

colors=colorRampPalette(brewer.pal(9,"Blues"))(7)


level1<-c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Lipidomics', 'Metagenomics', 'miRNAs')
level2<-c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Metagenomics')

# Process; if methylation or histone; add epigenomics!
preprocessing<-function(df,colname){
  #' Split the column 
  #' Return the split variables to get the frequencies 
  
  
  splitted<-str_split(df[[colname]], ',|\r|\n') # split by space, comma, newline
  splitted<-lapply(splitted,trimws)           # remove whitespace
  splitted<-splitted[!is.na(splitted)]      # remove nas 
  
  
  return(splitted)
  
}




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
  #'These groups are for objective - method
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*','.*understand molecular.*'),
                                  c('Diagnosis/Prognosis', 'understand molecular mechanisms')))
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





expand_ObjeMeth<- function(stats){
  # expand the objective method column!! 
  # this will increase the rows of stats dataframe and fill in the objective and method columns 
  stats<-stats %>% 
    mutate(ObjeMeth=strsplit(ObjeMeth, ',|\r|\n' ))%>%
    unnest(ObjeMeth) 
  
  stats <-stats %>% separate(ObjeMeth, c("objective","method"), sep = " - ")
  return(stats)
}

library('readxl')
library('stringr')
library(ggplot2)
library(data.table)
sysinf <- Sys.info()


# stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
# stats<-read_excel('H:/My Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )
os <- sysinf['sysname']
if ( os  == 'Darwin'){
  stats<-read_excel('/Users/efiathieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )
}else{
  stats<-read_excel('E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )
  }

###Filters
#### 1. remove same sample 
#stats_filter<-stats[stats$`Data` %like% 'Proteomics',]


#### Split the omics and count number each used 
#' 1. Split to same sample yes or no
#' 2 Simplify to lower level 
#' 3. Plot 
#' 
#' 

#stats$Cancer<-c(rep('no',345), rep('yes',(nrow(stats)-345)))
#stats=stats[1:600,] do not filter anymore 
### Remove reviews, remove rejected articles
## GLOBAL FILTER
stats <- stats %>%
  filter(Type!= 'Review' | is.na(Type)) %>%
  filter(is.na(`Rejection /Critic`)) %>%
  filter(tolower(same_sample)!='no' | is.na(same_sample)) %>% # also includes nas that i did not label as no
  filter(!is.na(`Objective-Method`))




stats$same_sample<-as.factor(tolower(stats$same_sample))
stats$Cancer<-as.factor(tolower(stats$Cancer))

x_group<-'Cancer'

get_frequencies_by_group<-function(stats,colname){
  
  df_by_group<- stats %>%
    group_by_at(x_group) %>%
    group_modify(~ preprocessing(.x, colname) %>%
                   get_frequencies() 
    )
  
  
  return(df_by_group)
}

stats$same_sample<-as.factor(tolower(stats$same_sample))



colname<-'Data'
#### Also filter by omics
frequencies_by_group<-get_frequencies_by_group(stats, colname)

freq_to_plot<-frequencies_by_group %>% filter(Var1 %in% tolower(level2))
single_omics_frequencies=freq_to_plot

#### Create the stacked plot

ggplot(stats[stats$Cancer=='no',])+stat_count(aes(tolower(Disease)))+
  theme(axis.text.x = element_text(size=rel(plt_txt_size),angle = 35, vjust = 0.5, hjust=1))
  
  ### TODO: how many are the total number considred!!??

library(grid)

plotByData<-function(df_by_group, y_group){
  ggplot(df_by_group, 
         aes(x=reorder(Var1, -Freq, sum), y=Freq))+
    aes_string(fill=factor(y_group)) +#, labels=c('Cancer', 'Other Diseases', 'NA')))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL, title=paste0('Combinations with > ',freq_cutoff, ' occurences'))+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 35, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,1.7,2.5),"cm"))
    
  ggsave(paste0('plots/byCombinations', as.character(colname), '.png'), width = 8, height=6)
  
}

plt_txt_size=1.5
plotByData(freq_to_plot, y_group=x_group)


  ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq))+
         aes_string(fill=x_group)+
  geom_bar(stat='identity',position='stack')+
  labs(x=NULL)+
  theme(axis.text.x = element_text(size=rel(plt_txt_size),angle = 25, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,2,1.5),"cm"))+
  labs(y='Frequency')+
   scale_fill_discrete(name = " ", labels = c("Other Diseases", "Cancer"))
    

ggsave(paste0('plots/SingleOmicsby', as.character(colname), '.png'), width = 6, height = 5)



####### Objectives 

##### Plot by objective

#### Get combinations 
### Co-Occurrence #
library(tidyverse)

colname<-'Data'

omics_data<-df_by_group[[1]]


####
#' Combinations and their objectives
#' 
#' 

get_combs<- function(x, omics_level=level1){
    
  #' return omics combinations as individual strings to count them+
  x<-unlist(x)
    x<-x[tolower(x) %in% tolower(omics_level)]
    if (length(x)>1){
      x<-x[order(x)]
      combn(x,2, FUN=paste, collapse=' - ')
    }
}

preprocessing_combinations(preprocessing(stats, 'Data'))
preprocessing_combinations(preprocessing(new, 'Data'))

preprocessing_combinations<-function(x){
  #' Create combinations of omics datasets  
  x<-x[x!='']
  x<-x[!is.na(x)]
  #' Create pairs of omics 
  #' #
  #'
  combinations<-lapply(x,get_combs, omics_level=level2)
  return(combinations)
}


total<-NROW(stats[!is.na(stats$Data),]$PMID)

y_group='same_sample'

y_group='Cancer'

new<-stats[! ( is.na(stats['Data'] ) ),]



df_by_group <- new %>%
  group_by_at(y_group) %>%
  group_modify(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations %>%
              get_frequencies() 
  )
freq_cutoff<-0
df_by_group_filtered<-df_by_group %>% 
  group_by(Var1)  %>% 
  filter( sum(Freq) >= freq_cutoff) 

p<-plotByData(df_by_group_filtered, y_group)

p

df_by_group$perc<-as.numeric(df_by_group$Freq)/(NROW(new))*100


cancer_filter=c('no')

df_by_group_data<-df_by_group 




plotGridCombinations<-function(df_by_group){
  
  #### Plots a grid showing the frequencies 
  combinations <-df_by_group %>% separate(Var1, c("Omics1","Omics2"), sep = " - ")
  combinations<-combinations[order(combinations$Freq, decreasing = TRUE),]
  
  
  g<-ggplot(combinations)+aes(Omics1, Omics2, fill=abs(Freq)) +
    geom_tile()+
    geom_text(aes(label = round(Freq, 2)), size=7)+
    theme(axis.text.x = element_text(size=rel(plt_txt_size),angle = 45, vjust = 0.5, hjust=1), 
                                     axis.text.y = element_text( size=rel(1.5)), 
          plot.margin = margin(10, 10, 40, 20))+
    labs(x=NULL, y=NULL)+
    scale_fill_gradient(low = "white", high = "red")+
    guides(fill=guide_legend(title="Frequency"))+
    facet_grid(~Cancer,  labeller = labeller(Cancer=
                                               c('no'='Other Diseases','yes' ='Cancer')),  
               scales = 'free', space='free')
  show(g)
  ggsave(paste0('plots/GridPlot', as.character(colname), '.png'), width = 10, height=6)
  
}

plotGridCombinations(df_by_group_data)


#TODO: add red find overflow 
  
  

comb_frequencies_by_group<-df_by_group



##### List of disease by combinations 



new_disease<-new
new_disease$disease_group<-group_disease(new_disease, 'Disease')$Disease
new_disease$disease_group[which(new_disease$Cancer=='yes')]
new_disease$disease_group[which(new_disease$Cancer=='yes')]<-'Cancer'

x_group='disease_group'

df_by_group<-new_disease %>%
  group_by_at(c(x_group, 'Cancer')) %>%
  group_modify(~ preprocessing(.x, colname)  %>%
                 preprocessing_combinations %>%
                 get_frequencies() 
  )  



df_by_group_disease<- df_by_group

#### do the filtering further down... jump to line 418
# TODO: make this automatic to create all graphs by a switch



show_p<-plotbyObjective(df_to_plot )

df_nested<-df_by_group[,-4] %>% 
  nest(data = disease_group)
# size of tibbles shows frequencies
df_nested$Freq<-sapply(df_nested$data, dim)[1,]

df_nested$concat<-sapply(df_nested$data, function(x){
  x2=as.character(unlist(x))
  print(as.character(x2))
  x2=x2[order(x2)]
  print(x2)
  x2=paste(x2, collapse=', ')
  return(x2)
})

df_nested %>% filter()
df_nested<-df_nested %>% arrange(Cancer, desc(Freq))
# two lists - most frequest is just not cancer now 
df_nested_filtered<- df_nested %>% filter(Var1 %in% most_frequent$Var1[1:7])
df_nested_filtered<-df_nested %>% filter(Cancer %in% c('yes', 'no'))%>%
                     group_by(Cancer) %>% 
                      slice_max(order_by = Freq, n=5)

write.table(df_nested[,-3], file = "review/output/data_diseases.txt", sep='\t', row.names = FALSE, quote = FALSE)






########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 
colnames(stats)[which(colnames(stats)=='Objective-Code')]<-'objective'
colnames(stats)[which(colnames(stats)=='Integration method-Category')]<-'method'
colnames(stats)[which(colnames(stats)=='Objective-Method')]<-'ObjeMeth'



stats_expanded<-expand_ObjeMeth(stats)


stats_fil<-stats_expanded[stats_expanded$Cancer == cancer_filter,]
stats_fil<-stats_expanded

#stats_fil$disease_group<-group_disease(stats_fil$Disease)

#change here to select by objective or by method 
# TODO: MAKE this one variable to choose method or objective






x_group<-'method'
# and here!!! 
new<-stats_fil %>% 
  group_by('Cancer') %>%
  
  mutate(method=strsplit(method, ',|\r|\n' ))%>%
  unnest(method) 


x_group<-'objective'
new<-stats_fil %>%
  group_by('Cancer') %>%
  mutate(objective=strsplit(objective, ',|\r|\n' ))%>%
  unnest(objective)


colname='Data'
new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))
new<-group_objectives(new, 'objective')






#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
# x group is the objective
df_by_group<-new %>%
  group_by_at(c(x_group, 'Cancer')) %>%
  group_modify(~ preprocessing(.x, colname)  %>%
              preprocessing_combinations %>%
              get_frequencies() 
  )  




# Attach the key names back to the dataframe 
df_by_group<-as.data.frame(as.matrix(df_by_group))
df_by_group['key_names']<-df_by_group[x_group]
df_by_group$Freq<-as.numeric(df_by_group$Freq)


# non cancer: 10,7, cancer: 7,3,

df_to_plot<-df_by_group



##TODO: MOVE TO FUNCTION
#overwrite


x_group<-'disease_group'
if (x_group == 'objective'){
  df_most_common<-filter_common_groups(df_by_group, freq_cutoff = c(25,25))
  df_to_plot<-df_most_common
  df_to_plot<-relabel_objectives_short(df_to_plot)
  df_to_plot<-df_to_plot[!df_to_plot[x_group]=='NA',]
  plot_width=10
  plot_height=9
  
}else if (x_group == 'disease_group' ){
  # select common groups 
  df_by_group<-df_by_group_disease
  df_most_common<-filter_common_groups(df_by_group,  freq_cutoff = c(9,6))
  
  df_to_plot<-df_most_common
  # show them all
  #df_to_plot<-df_by_group
  
  # remove cancer
  df_to_plot<-df_to_plot[!df_to_plot[x_group]=='NA',]
  df_to_plot<-df_to_plot[!df_to_plot[x_group]=='all diseases',]
  
  df_to_plot<-df_to_plot[!df_to_plot['Cancer']=='NA',]
  
  plot_width=9
  plot_height=5.5
  
}


# Remove non_cancer
df_to_plot=  plot_filters(df_to_plot)

show_p<-plotbyObjective(df_to_plot, plot_width=plot_width, plot_height = plot_height, plot_cols = TRUE)
show_p







##### 
#' Create an edge list from the frequency table ! 
# aggregated or separately? 

## todo label the size of the nodes df_to_plotby frequency 
comb_freq<- comb_frequencies_by_group %>% filter(Cancer ==cancer_filter)
single_omics_frequencies_filtered<-single_omics_frequencies %>% filter(Cancer ==cancer_filter)

edge_list<-data.frame(do.call(rbind, str_split(comb_freq$Var1, ' - ')))
edge_list$weight<-comb_freq$Freq
edge_list<-edge_list[order(edge_list$weight, decreasing = TRUE),]

# most frequent datasets
# most_frequent<-comb_freq[order(comb_freq$perc, 
#                                                decreasing = TRUE),]
# most_frequent_5<-most_frequent$Var1[1:5]
# edge_list

library(igraph)

#### Add the frequency of single omics to the network vertices
df<-single_omics_frequencies_filtered

net<-graph_from_data_frame(edge_list, directed = FALSE, vertices =NULL)
net_att<-df[match(V(net)$name, df$Var1),]
cancer_filter

# TODO: assign omics frequencies of all samples or only cancer? 
#vertex_attr(net, 'freq', index=V(net))<-single_omics_frequencies$Freq

p<-plot.igraph(net, edge.width=edge_list$weight, vertex.size=log2(net_att$Freq)*5)
#save(p,paste0('plots/network', as.character(colname), '.png'))


#g <- set.vertex.attribute(g,'id',1,'first_id')












