# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics

colname<-'Data'
library('dplyr')
library('purrr')
level1<-c('transcriptomics', 'genomics','epigenomics', 'proteomics', 'metabolomics', 'metagenomics', 'mirnas')

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



group_objectives_method<-function(df, Var1){
  #'Group objective code column 
  #'These groups are for objective - method
  df[Var1]<-sapply(df[Var1],
                   function(x) 
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*','.*understand mol.*'),
                           c('Diagnosis/Prognosis', 'understand molecular mechanisms')))
  return(df)
}


df<-new
Var1<-'method'
library(gsubfn)
group_methods<-function(df, Var1){
      df[Var1]<-sapply(df[Var1],function(x){
                 mgsub::mgsub(tolower(x),  
                c(".*learning.*|.*decision.*|.*neural.*|.*boosting.*|.*kmeans.*|.*support vector.*",  
                  '.*pca.*|.*cluster.*', '.*regression.*|.*linear model.*', '.*factor.*|.*decomposition.*', 
                  '.*multivar.*', '.*snf.*|.*network.*', '.*gsea.*', '.*cca.*', 
                  '.*kernel.*', '.*autoencoder.*', 
                  '.*partial least.*|.*diablo.*'), 
                c( "machine/deep learning", 'clustering',
                                 'regression', 'factor analysis', 'multivariate analysis', 
                   'network', 'enrichment', 'canonical correlation analysis',
                   'kernel learning', 'autoencoder', 'partial least squares'
                   ))}
)
      #new_col=as.factor(new_col)
return(df)                 
}
library(magrittr)



library('readxl')
library('stringr')
library(ggplot2)
library(data.table)



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
stats<-read_excel('E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )
stats<-read_excel('/Users/efiathieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )
#stats<-stats[1:289,]
stats$PMID<-as.numeric(stats$PMID)
stats$Cancer<-c(rep('no',345), rep('yes',(nrow(stats)-345)))

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

cancer_filter = 'no'
new<-new[new$Cancer == cancer_filter,]


new <-new %>% separate(ObjeMeth, c("objective","method"), sep = " - ")

# this can produce the plot of x axis objective with groups of method 
# or the other way round! 

#x_group<-'method'
#colname='objective'
width=10

x_group<-'objective'
colname='method'
width=7



new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))
new[colname] <-apply(new[colname], 1, function(x) trimws(tolower(x)))


new[x_group]

new<-group_objectives_method(new, 'objective')
new<-group_methods(new, 'method')


keys<-pull(new %>%
             group_by_at(x_group) %>%
             group_keys())




#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
new<-new[! (is.na(new[x_group]) | is.na(new['Data'] )),]
df_by_group<-new %>%
  group_by_at(x_group) 

df_by_group<-new %>%
  group_by_at(x_group) %>%
  group_modify(~ preprocessing(.x, colname) %>%
              get_frequencies() 
  )  
# %>%  map_df(I, .id=x_group) 



# Attach the key names back to the dataframe 
df_by_group<-as.data.frame(as.matrix(df_by_group))

df_by_group['key_names']<-df_by_group[x_group]

df_by_group$Freq<-as.numeric(df_by_group$Freq)
df_by_group$perc<-as.numeric(df_by_group$Freq)/(NROW(new))*100

df_to_plot<-df_by_group %>%
  group_by(Var1)  %>%
  filter( sum(Freq) >= 3) %>%
  group_by_at(x_group)  %>%
  filter( sum(Freq) >= 2)


show_p<-plotbyObjective(df_to_plot )


show_p

plotbyObjectie<-function(df){ 
  g<-ggplot(df, aes(x=reorder(key_names, -Freq, sum), y=perc, fill=Var1))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL)+
    theme(axis.text.x = element_text(size=rel(1.5),angle = 25, vjust = 0.5, hjust=1))+

    theme(plot.margin=unit(c(1,1,3,3),"cm"))
  ggsave(paste0('plots/byObjMethod', as.character(x_group),'_', cancer_filter, '.png'), width = width, height=6)
  return(g)
  
  
}

new_concise<-new_concise[!is.na(new_concise['Data']),]
new_concise<-new[c('Data', 'objective', 'method' )]




#install.packages('alluvial')
#install.packages('ggalluvial')
install.packages('ggsankey')
library('ggalluvial')
library('alluvial')

library('ggsankey')


new2<-new %>% 
  mutate(Data=strsplit(Data, ',|\r|\n' ) )%>%
  unnest(Data) 

cancer_filter=c("yes")
new2<-new2 %>% filter(Cancer %in% cancer_filter)


new2$Data<-tolower(trimws(new2$Data))
new2<-new2 %>% filter(Data %in% tolower(level1))

levels(as.factor(new2$Data))

new2<-new2[!is.na(new2$method),]
new2<-new2[!is.na(new2$Data),]
new2<-new2[!is.na(new2$objective),]

axis1='objective'
axis2='method'

counts<-new2 %>% count(objective, method)


counts<-counts%>% filter(n>1)


df<-counts
ggplot(as.data.frame(df),
       aes_string(y = 'n', axis1 = axis1, axis2 = axis2)) +
  geom_alluvium(aes_string(fill = axis1),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) + 
  
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:2, labels = c(axis1,axis2)) +
  ggtitle(paste0("Multi omics objectives, Cancer = ", cancer_filter))


ggsave(paste0('plots/alluvial', as.character(paste0(axis1, axis2)),'_', cancer_filter, '.png'), width = 7, height=6)



counts<-new2 %>% count(objective, method)
counts<-counts%>% filter(n>1)




# New implementation with ggalluvial
axis1='objective'
axis2="method"
counts <- new2 %>% 
  count( objective, method) %>% 
  mutate(
    col = objective
  ) %>%
  ggalluvial::to_lodes_form(key = type, axes = c(axis1, axis2))

df<-counts %>% filter(n>2)
# df<-counts
ggplot(data = df, aes(x = type, stratum = stratum, alluvium = alluvium, y = n)) +
  # geom_lode(width = 1/6) +
  geom_flow(aes(fill = col), width = 1/6, color = "darkgray",
            curve_type = "cubic") +
  # geom_alluvium(aes(fill = stratum)) +
  geom_stratum(color = "grey", width = 1/6) + 
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  theme(
    panel.background = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15, face = "bold"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_viridis_d()+
  ggtitle(paste0("Multi omics objectives, Cancer = ", cancer_filter))


ggsave(paste0('plots/ggalluvial', as.character(paste0(axis1, axis2)),'_', cancer_filter, '.png'), width = 7, height=6)



