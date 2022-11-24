# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics
colname<-'Data'
library('dplyr')
library('purrr')
library(RColorBrewer)
source('utils.R')
library("viridis")           # Load



### Load the package or install if not present
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

colors=colorRampPalette(brewer.pal(9,"Blues"))(7)


level1<-c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Lipidomics', 'Metagenomics', 'miRNAs')
level2<-c('Transcriptomics', 'Genomics','Epigenomics', 'Proteomics', 'Metabolomics', 'Metagenomics')


remove_objectives<-c('multiomics pathway analysis', 'biomarker discovery', 'other', 
                     'molecular interactions', 'understand molecular mechanisms', 
                     'downstream', 
                     'drug repurposing', 'missing data')


# Process; if methylation or histone; add epigenomics!
preprocessing<-function(df,colname){
  #' Split the column 
  #' Return the split variables to get the frequencies 
  
  
  splitted<-str_split(df[[colname]], ',|\r|\n') # split by space, comma, newline
  splitted<-lapply(splitted,trimws)           # remove whitespace
  splitted<-splitted[!is.na(splitted)]      # remove nas 
  
  
  return(splitted)
  
}


freq_cutoff<-0

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
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*',
                                               '.*connect.*',
                                               'multiomics pathway analysis'),
                                  c('diagnosis/ prognosis', 
                                    'Extract complex patterns', 
                                    'downstream')))
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
  stats<-read_excel('/Users/efiathieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge_2.xlsx' )
}else{
  stats<-read_excel('D:/DATADRIVE/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge_2.xlsx' )
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

#' remove nas - careful check and fill in if you forgot to add cancer status 
stats[is.na(stats$Cancer),'PMID']



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
  theme(axis.text.x = element_text(size=rel(rel_txt),angle = 35, vjust = 0.5, hjust=1))
  
  ### TODO: how many are the total number considred!!??

library(grid)

plotByData<-function(df_by_group, y_group){
  ggplot(df_by_group, 
         aes(x=reorder(Var1, -Freq, sum), y=Freq))+
    aes_string(fill=factor(y_group)) +#, labels=c('Cancer', 'Other Diseases', 'NA')))+
    geom_bar(stat='identity',position='stack')+
    labs(x=NULL, title=paste0('Combinations with > ',freq_cutoff, ' occurences'))+
    theme(axis.text.x = element_text(size=rel(1.3),angle = 35, vjust = 0.5, hjust=1))+
    theme(plot.margin=unit(c(1,1,1.7,2.5),"cm") )+
    scale_color_viridis(option = "D")
    
  ggsave(paste0('plots/byCombinations', as.character(colname), '.png'), width = 8, height=6)
  
}

rel_txt=1.5
plotByData(freq_to_plot, y_group=x_group)


freq_to_plot<-freq_to_plot[!is.na(freq_to_plot$Cancer),]
pal_npg("nrc", alpha = 0.9)(10)

freq_to_plot<-freq_to_plot[!is.na(freq_to_plot$Cancer),]
  p<- ggplot(freq_to_plot, aes(x=reorder(Var1, -Freq, sum), y=Freq))+
         aes_string(fill=x_group)+
  geom_bar(stat='identity',position='dodge',  color='black')+
  labs(x=NULL)+
  theme(axis.text.x = element_text(size=rel(rel_txt),angle = 25, vjust = 0.5, hjust=1))+
  theme(plot.margin=unit(c(1,1,2,1.5),"cm"))+
  labs(y='Frequency')+
   scale_fill_discrete(name = " ", labels = c("Other Diseases", "Cancer"))+
    scale_fill_manual(values =c('#F39B7FE5', '#8491B4E5'))
    
show(p)
  ggsave(paste0('plots/SingleOmicsby', as.character(colname), '.png'), width = 6, height = 5)



####### Objectives 

##### Plot by objective

#### Get combinations 
### Co-Occurrence #
library(tidyverse)

colname<-'Data'



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
preprocessing_combinations(preprocessing(stats, 'Data'))
preprocessing_combinations(preprocessing(new, 'Data'))

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

df_by_group_filtered<-df_by_group %>% 
  group_by(Var1)  %>% 
  filter( sum(Freq) >= freq_cutoff) 

p<-plotByData(df_by_group_filtered, y_group)

p

df_by_group$perc<-as.numeric(df_by_group$Freq)/(NROW(new))*100


cancer_filter=c('no')

df_by_group_data<-df_by_group 

text_size=16


plotGridCombinations<-function(df_by_group){
  
  #### Plots a grid showing the frequencies 
  combinations <-df_by_group %>% separate(Var1, c("Omics1","Omics2"), sep = " - ")
  combinations<-combinations[order(combinations$Freq, decreasing = TRUE),]
  
  
  g<-ggplot(combinations)+aes(Omics1, Omics2, fill=abs(Freq)) +
    geom_tile()+
    geom_text(aes(label = round(Freq, 2)), size=rel(6))+
    theme(axis.text.x = element_text(size=text_size,angle = 45, vjust = 0.5, hjust=1), 
                                     axis.text.y = element_text( size=text_size), 
          plot.margin = margin(10, 10, 40, 20))+
    labs(x=NULL, y=NULL)+
    scale_fill_gradient(low = "white", high = "red")+
    guides(fill=guide_legend(title="Frequency"))+
    facet_grid(~Cancer,  labeller = labeller(Cancer=
                                               c('no'='Other Diseases','yes' ='Cancer'),size=text_size),  
               scales = 'free', space='free')+
    theme(strip.text.x = element_text(size = text_size))
  show(g)
  ggsave(paste0('plots/GridPlot', as.character(colname), '.jpeg'), width = 10, height=6)
  
}
df_by_group_data<-df_by_group_data[!is.na(df_by_group_data$Cancer),]
plotGridCombinations(df_by_group_data)


#TODO: add red find overflow 
  
  

comb_frequencies_by_group<-df_by_group_data



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
#df_nested_filtered<- df_nested %>% filter(Var1 %in% most_frequent$Var1[1:7])
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




new[c('objective')]

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
df_by_group$Freq<-as.numeric(df_by_group$Freq)




    

add_percentage<-function(df){
  
  df<-df %>%
    group_by_at(c(x_group, 'Cancer')) %>%
    mutate(percent=Freq/sum(Freq)*100)
  
  return(df)
}

##TODO: MOVE TO FUNCTION
#overwrite

# Switch here for both 
x_group<-'objective'
#x_group<-'disease_group'
top_n<-c(8,8)

 if (x_group=='disease_group'){
  df_by_group<-df_by_group_disease
  top_n<-c(8,11)
  
  
}
  
  


###### Group the not common combinations together and label as none
most_common<-filter_common_groups(df_by_group, top_n =top_n, x_group=x_group)
most_common_pairs<-unlist(most_common[[1]])
most_common_groups<-unlist(most_common[[2]])
most_common_pairs
most_common_groups; length(most_common_groups)
freq_cutoff_objectives<-c(5,5)
df_filt<-df_by_group
df_filt$Var1<-as.factor(df_filt$Var1)

levels(df_filt$Var1)[!(levels(df_filt$Var1) %in% most_common_pairs)]<-'other'
## Regroup and sum the frequency  combinations after renaming
df_filt<-df_filt %>% 
  group_by_at(c(x_group, 'Cancer', 'Var1')) %>%
  summarise(Freq=sum(Freq), .groups = 'drop') %>% 
  as.data.frame() 

# put the other category last! 
new_l<-levels(df_filt$Var1);
new_l<-c(  'other',new_l[-which(new_l=='other')])
df_filt$Var1<-factor(df_filt$Var1, levels=new_l)
df_filt['key_names']<-df_filt[x_group]

# Remove low frequency groups
df_filt<-df_filt[(df_filt[,x_group] %in% c(most_common_groups)),]

if (x_group == 'objective'){
  
  df_to_plot<-df_filt
  df_to_plot<-relabel_objectives_short(df_to_plot)
  df_to_plot$Var1<-relabel_omics_short(df_to_plot$Var1)
  
  df_to_plot<-df_to_plot[!df_to_plot[x_group]=='NA',]

  df_to_plot<-df_to_plot[!(df_to_plot$objective %in% remove_objectives),]
  
  df_to_plot<- add_percentage(df_to_plot)
  
  plot_width=12
  plot_height=3.5
  plot_cols=TRUE
  stack_horiz=TRUE
  xangle=0
  
  }else if (x_group == 'disease_group' ){
  # select common groups 

  df_to_plot<-df_filt
  # show them all
  #df_to_plot<-df_by_group
  
  df_to_plot$Var1<-relabel_omics_short(df_to_plot$Var1)
  
  # remove cancer
  df_to_plot<-df_to_plot[!df_to_plot[x_group]=='NA',]
  df_to_plot<-df_to_plot[!df_to_plot[x_group]=='all diseases',]
  
  df_to_plot<-df_to_plot[!df_to_plot['Cancer']=='NA',]
  df_to_plot<- add_percentage(df_to_plot)
  
  
  
  plot_width=8
  plot_height=3.5
  plot_cols=FALSE
  stack_horiz=TRUE

  
}


# Remove non_cancer

df_to_plot<-df_to_plot[df_to_plot$Cancer %in% c('yes', 'no'),]
# df_to_plot=  plot_filters(df_to_plot)




show_p<-plotbyObjective(df_to_plot, plot_width=plot_width, plot_height = plot_height, plot_cols = plot_cols,
                        stack_horiz = stack_horiz)
show_p

# df_to_plot[df_to_plot$disease_group=='Metabolism',]

length(most_common_pairs)








