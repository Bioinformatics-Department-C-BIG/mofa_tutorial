# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics

#install.packages('purrr')

colname<-'Data'
library('dplyr')
library('purrr')
source('literature_review.R')
source('utils.R')


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
Var1<-'method'
#install.packages('mgsub')
library(gsubfn)
library(mgsub)


library(magrittr)

#install.packages('readxl')
library('readxl')
library('stringr')
library(ggplot2)
#install.packages('data.table')
library(data.table)

#install.packages('gdata')
library(gdata) 



sysinf <- Sys.info()


# stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
# stats<-read_excel('H:/My Drive/PHD 2020/Literature/Data Integration/Multi-omics_not cancer_merge.xlsx' )
os <- sysinf['sysname']
if ( os  == 'Darwin'){
  data_int_dir<-'/Users/efiathieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/'
  }else{
  data_int_dir<-'E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/'
  
}


stats<-read_excel(paste0(data_int_dir,'Multi-omics_merge.xlsx' ))




#stats<-stats[1:600,]
stats$PMID<-as.numeric(stats$PMID)
#stats$Cancer<-c(rep('no',345), rep('yes',(nrow(stats)-345)))


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
  filter(Type!= 'Review' | is.na(Type)) %>%
  filter(is.na(`Rejection /Critic`)) %>%
  filter(tolower(same_sample)!='no' | is.na(same_sample)) %>%
  filter(!is.na(Data))

#### Create the stacked plot


#### Create the stacked plot

library(grid)
#install.packages('tidyverse')
library(tidyverse)

 ########
###
#' Expand the objective-code
#' And get frequencies by objective group 
#' 
colnames(stats)[which(colnames(stats)=='Objective-Code')]<-'objective'
colnames(stats)[which(colnames(stats)=='Integration method-Category')]<-'method'
colnames(stats)[which(colnames(stats)=='Objective-Method')]<-'ObjeMeth'


cancer_filter = 'yes'


new<-expand_ObjeMeth(stats)


width=10

x_group<-'objective'
colname='method'
width=7


new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))
new[colname] <-apply(new[colname], 1, function(x) trimws(tolower(x)))


new<-group_objectives_method(new, 'objective')
# make a backup of the original entries 
new['method_orig']<-new['method']
new<-group_methods(new, 'method')

as.data.frame(new[new$method=='regression',c('PMID', 'method_orig')])
#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
new<-new[! (is.na(new[x_group]) | is.na(new['Data'] )),]


new_objective_method<-new

df_by_group<-new %>%
  group_by_at(c(x_group, 'Cancer')) %>%
  group_modify(~ preprocessing(.x, colname) %>%
              get_frequencies() 
  )  


# Attach the key names back to the dataframe 
df_by_group<-as.data.frame(as.matrix(df_by_group))



df_by_group$Freq<-as.numeric(df_by_group$Freq)
df_by_group$perc<-as.numeric(df_by_group$Freq)/(NROW(new))*100




######## Plotting - Filter
df_most_common<-filter_common_groups(df_by_group,freq_cutoff=c(5,5))
# group all the miscellaneous in one category
most_common_groups<-levels(as.factor(df_most_common$Var1))
most_common_objectives<-levels(as.factor(df_most_common$objective))
df_to_plot<- df_by_group

other_tools<-df_to_plot$Var1[!(df_to_plot$Var1 %in% most_common_groups)]; 
other_tools
df_to_plot$Var1[!(df_to_plot$Var1 %in% most_common_groups)]<-'Other'
other_objectives<-df_to_plot$objective[!(df_to_plot$objective %in% most_common_objectives)]
other_objectives


df_to_plot$objective[!(df_to_plot$objective %in% most_common_objectives)]<-'Other biomarker discovery questions'

# remove most common objective -only for visualization
df_to_plot<-df_to_plot[(df_to_plot$objective %in% most_common_objectives),]
df_to_plot<-df_to_plot[!(df_to_plot$objective == 'multiomics pathway analysis'),]

df_to_plot<-df_to_plot[(df_to_plot$Var1 %in% most_common_groups),]
df_to_plot<-df_to_plot[!(df_to_plot$Var1 == 'multiomics pathway analysis'),]


df_to_plot$Var1<-as.factor(df_to_plot$Var1)
df_to_plot=df_to_plot[df_to_plot$Cancer %in% c('yes', 'no'),]
df_to_plot['key_names']<-df_to_plot[x_group]

df_to_plot<-relabel_objectives_short(df_to_plot)

#df_to_plot %>% 
 # group_by(c(Cancer, objective, key_names, Var1)) 

show_p<-plotbyObjective(df_to_plot, 'Methods', plot_width=9, plot_height=9)

# TODO REGROUP THE TOTALS 

show_p

new_concise<-new_concise[!is.na(new_concise['Data']),]
new_concise<-new[c('Data', 'objective', 'method' )]





















axis1='objective'
axis2="Data"

df_to_plot[axis2]<-df_to_plot$Var1
df_to_plot<-df_to_plot%>%filter(Cancer==cancer_filter)
df_to_plot$n=df_to_plot$Freq
df_to_plot$col<-df_to_plot$key_names
run_sankey(df_to_plot, axis1,axis2, cancer_filter  )







################ Section 5: ALLUVIAL 


#install.packages('alluvial')
#install.packages('ggalluvial')
#install.packages('ggsankey')
library('ggalluvial')
library('alluvial')

library('ggsankey')

cancer_filter=c("yes")
new2<-new %>% 
  mutate(Data=strsplit(Data, ',|\r|\n' ) )%>%
  unnest(Data) 


new2<-new %>% filter(Cancer %in% cancer_filter)


new2$Data<-tolower(trimws(new2$Data))
new2<-new2 %>% filter(Data %in% tolower(level1))

levels(as.factor(new2$Data))


new2<-new2[!is.na(new2$method),]
new2<-new2[!is.na(new2$Data),]

new2<-new2[!is.na(new2$objective),]

axis1='Data'
axis2='objective'


counts<-new2 %>% count(Data, objective)

if (cancer_filter == 'yes')
  {n_cutoff=1} else {n_cutoff=1}

# New implementation with ggalluvial
axis1='Data'
axis2='objective'

counts1 <- new2 %>% 
  count( Data, objective) %>% 
  mutate(
    col = objective
  ) %>%
  ggalluvial::to_lodes_form(key = type, axes = c(axis1, axis2))




run_sankey(new2, axis1, axis2, cancer_filter  )



  stats
  
  
  # Extract information 
  ## TODO: CREATE a new column whch is merged from objective and method
  # TODO: add a \newline after each objective 
    stats_to_write<-stats[c('PMID','Title', 'Data', 'ObjeMeth')]
  stats_to_write$x<-'\\\\'
  
  stats_to_write$ObjeMeth<- gsub(',',',\\\\newline',stats_to_write$ObjeMeth)
  #install.packages('gdata')

write.table(stats_to_write, file = "review/output/literature_latex_2.txt", sep = " & ", row.names = FALSE, quote = FALSE)



