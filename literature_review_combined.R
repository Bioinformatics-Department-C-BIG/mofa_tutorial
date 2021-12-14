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
                     gsub('.*diagnosis.*|*prognosis*', 'Diagnosis/Prognosis', tolower(x)))
  return(df)
}



df<-new
Var1<-'method'
library(gsubfn)
group_methods<-function(df, Var1){
      new_col<-sapply(df[Var1],function(x){
                 mgsub::mgsub(tolower(x),  
                c(".*learning.*|.*decision.*|.*neural.*",  '.*pca.*', '.*regression.*', '.*factor.*', 
                  '.*multivar.*', '.*snf.*', '.*gsea.*', '.*cca.*', 
                  '.*kernel.*'), 
                c( "machine/deep learning", 'clustering',
                                 'regression', 'factor Analysis', 'multivariate analysis', 
                   'network', 'enrichment', 'canonical correlation analysis',
                   'kernel learning'
                   ))}
)
      #new_col=as.factor(new_col)
return(new_col)                 
}
library(magrittr)



library('readxl')
library('stringr')
library(ggplot2)
library(data.table)



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
stats<-read_excel('G:/My Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )
stats$PMID<-as.numeric(stats$PMID)

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
  #filter(Type!= 'Review')%>%
  filter(is.na(`Rejection /Critic`))

which(stats$PMID==31856727)

#### Create the stacked plot

library(grid)
library(tidyverse)

########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 
colnames(stats)[which(colnames(stats)=='Objective-Code')]<-'objective'
colnames(stats)[which(colnames(stats)=='Integration method-Category')]<-'method'



new<-stats %>% 
  mutate(method=strsplit(method, ',|\r|\n' ))%>%
  unnest(method) 

#x_group<-'objective'
x_group<-'method'

new<-new %>%
 mutate(objective=strsplit(objective, ',|\r|\n' ))%>%
 unnest(objective)

colname='objective'

new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))


new<-group_objectives(new, 'objective')

new<-group_methods(new, 'method')
#which(new$PMID==31856727)


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
  group_map(~ preprocessing(.x, colname) %>%
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
  filter( sum(Freq) >= 2) %>%
  group_by_at(x_group)  %>%
  filter( sum(Freq) >= 1)

#df_to_plot<-df_by_group
df_to_plot<-df_to_plot[!is.na(df_to_plot$key_names),]

show_p<-plotbyObjective(df_to_plot )


show_p

plotbyObjective<-function(df){ 
  g<-ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=Freq, fill=Var1))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,2,1.7),"cm"))
  
  ggsave(paste0('plots/byGroupGroup', as.character(x_group), '.png'), width = 8, height=6)
  return(g)
  
  
}




#### Create the stacked plot

library(grid)
library(tidyverse)

########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 
colnames(stats)[which(colnames(stats)=='Objective-Code')]<-'objective'
colnames(stats)[which(colnames(stats)=='Integration method-Category')]<-'method'
colnames(stats)[which(colnames(stats)=='Objective-Method')]<-'ObjeMeth'



new<-stats %>% 
  mutate(ObjeMeth=strsplit(ObjeMeth, ',|\r|\n' ))%>%
  unnest(ObjeMeth) 



new <-new %>% separate(ObjeMeth, c("objective","method"), sep = " - ")



#x_group<-'objective'
x_group<-'method'
colname='objective'

new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))


new[x_group]

new<-group_objectives(new, 'objective')
new['method']<-group_methods(new, 'method')

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
  group_map(~ preprocessing(.x, colname) %>%
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
  filter( sum(Freq) >= 1) %>%
  group_by_at(x_group)  %>%
  filter( sum(Freq) >= 1)

#df_to_plot<-df_by_group
df_to_plot<-df_to_plot[!is.na(df_to_plot$key_names),]

show_p<-plotbyObjective(df_to_plot )


show_p

plotbyObjective<-function(df){ 
  g<-ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=Freq, fill=Var1))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 25, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,2,1.7),"cm"))
  
  ggsave(paste0('plots/byObjMethod', as.character(x_group), '.png'), width = 10, height=6)
  return(g)
  
  
}

new_concise<-new_concise[!is.na(new_concise['Data']),]
new_concise<-new[c('Data', 'objective', 'method' )]




#install.packages('alluvial')
#install.packages('ggsankey')




new_concise<-new[c('Data', 'objective', 'method', 'PMID' )]

new_concise<-new_concise %>% 
  mutate(Data=strsplit(Data, ',|\r|\n' ) )%>%
  mutate(Data=trimws(Data)) %>%
  mutate(Data=tolower(Data)) %>%
  unnest(Data) 

gather(new_concise) %>% 
  count(PMID) %>% 
  group_by(Data, method, objective)


alluvial

df<-df_to_plot
ggplot(as.data.frame(df),
       aes(y = Freq,
           axis1 = Survived, axis2 = Sex, axis3 = Class)) +
  geom_alluvium(aes(fill = Class),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:3, labels = c("Survived", "Sex", "Class")) +
  coord_flip() +
  ggtitle("Titanic survival by class and sex")






