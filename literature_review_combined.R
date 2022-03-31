# import literature review table and analysis 
# frequency of each disease 
# frequency of each omics

#install.packages('purrr')

colname<-'Data'
library('dplyr')
library('purrr')
source('literature_review.R')
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
                     mgsub::mgsub(tolower(x),c('.*diagnosis.*|*prognosis*','.*understand molecular.*'),
                           c('Diagnosis/Prognosis', 'understand molecular mechanisms')))
  return(df)
}


df<-new
Var1<-'method'
#install.packages('mgsub')
library(gsubfn)
library(mgsub)

group_methods<-function(df, Var1){
      df[Var1]<-sapply(df[Var1],function(x){
                 mgsub::mgsub(tolower(x),  
                c(".*learning.*|.*decision.*|.*neural.*|.*deep.*|.*autoencoder.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*cnn.*|.*classifier.*",  
                  '.*pca.*|.*cluster.*|.*kmeans.*|.*pins.*|.*movis.*|.*k means.*|.*partitioning.*', 
                  # '|.*lasso.*', 
                  '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                  '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*', 
                   '.*snf.*|.*net.*|.*piumet.*|.*omics integrator.*',
                  '.*gsea.*|.*enrichment.*', 
                  '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                  '.*correlation.*', 
                  '.*kernel.*', 
                  '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                  '.*ipa.*|.*activepathways.*|.*pathwaypca.*',
                  '.*david.*|.*metaboanalyst.*|.*nemo.*|.*adas.*'), 
                c( "ML/DL Classification", 
                   'ML/DL clustering',
                   'regression', 
                   'factor analysis', 
                   'network', 
                   'enrichment',
                   'correlation-based',
                   'correlation',
                   'kernel learning', 
                   'partial least squares',
                   'multiomics pathway analysis',
                   'Other tools'
                   ))}
)
      #new_col=as.factor(new_col)
return(df)                 
}



group_methods<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*neural.*|.*deep.*|.*autoencoder.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*cnn.*|.*classifier.*",  
                   '.*pca.*|.*cluster.*|.*kmeans.*|.*pins.*|.*movis.*|.*k means.*|.*partitioning.*', 
                   # '|.*lasso.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*', 
                   '.*snf.*|.*net.*|.*piumet.*|.*omics integrator.*',
                   '.*gsea.*|.*enrichment.*', 
                   '.*cca.*|.*smccnet.*|.*canonical correlation analysis.*',
                   '.*correlation.*', 
                   '.*kernel.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*',
                   '.*david.*|.*metaboanalyst.*|.*nemo.*|.*adas.*'), 
                 c( "ML/DL Classification", 
                    'ML/DL clustering',
                    'regression', 
                    'factor analysis', 
                    'network', 
                    'enrichment',
                    'correlation-based',
                    'correlation',
                    'kernel learning', 
                    'partial least squares',
                    'multiomics pathway analysis',
                    'Other tools'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}


group_methods<-function(df, Var1){
  df[Var1]<-sapply(df[Var1],function(x){
    mgsub::mgsub(tolower(x),  
                 c(".*learning.*|.*decision.*|.*boost.*|.*support vector.*|.*svm.*|.*random forest.*|.*cnn.*|.*classifier.*",  
                   ".*autoencoder.*|.*umap.*|.*tsne.*|.*deep.*.|*neural.*",
                   '.*pca.*|.*cluster.*|.*kmeans.*|.*pins.*|.*movis.*|.*k means.*|.*partitioning.*', 
                   # '|.*lasso.*', 
                   '.*regression.*|.*linear model.*|.*multivar.*|.*lasso.*|.*elastic net.*',
                   '.*factor.*|.*decomposition.*|.*mofa.*|.*intnmf.*', 
                   '.*snf.*|.*net.*|.*piumet.*|.*omics integrator.*',
                   '.*gsea.*|.*enrichment.*', 
                   '.*cca.*|.*smccnet.*|.*canonical correlation.*',
                   '.*correlation.*', 
                   '.*kernel.*', 
                   '.*partial least.*|.*diablo.*|.*pls.*|.*pls-da.*', 
                   '.*ipa.*|.*activepathways.*|.*pathwaypca.*',
                   '.*david.*|.*metaboanalyst.*|.*nemo.*|.*adas.*'), 
                 c( "ML Classification", 
                    'JDR - NonLinear',
                    'JDR - Linear',
                    'Regression', 
                    'JDR - Linear', 
                    'Graph-based', 
                    'enrichment',
                    'JDR - correlation-based',
                    'correlation',
                    'kernel learning', 
                    'JDR - partial least squares',
                    'multiomics pathway analysis',
                    'Other tools'
                 ))}
  )
  #new_col=as.factor(new_col)
  return(df)                 
}


library(magrittr)

#install.packages('readxl')
library('readxl')
library('stringr')
library(ggplot2)
#install.packages('data.table')
library(data.table)

#install.packages('gdata')
library(gdata) 



stats<-read_excel('C:/Users/athienitie/Google Drive/PHD 2020/Literature/Data Integration/Copy of Multi-omics_not cancer_updated at home  - November 2, 6_24 Pm.xlsx' )
stats<-read_excel('E:/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )
stats<-read_excel('/Users/efiathieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/Multi-omics_merge.xlsx' )
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
# new<-stats %>% 
#   mutate(ObjeMeth=strsplit(ObjeMeth, ',|\r|\n' ))%>%
#   unnest(ObjeMeth) 
# 
# new<-new[new$Cancer == cancer_filter,]
# 
# 
# new <-new %>% separate(ObjeMeth, c("objective","method"), sep = " - ")

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
df_most_common<-df_by_group %>%
  group_by(Var1, Cancer)  %>%
  filter( sum(Freq) >=3) %>%
  group_by_at(x_group)  %>%
  filter( sum(Freq) >= 5)

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
df_to_plot<-df_to_plot[(df_to_plot$Var1 %in% most_common_groups),]


df_to_plot$Var1<-as.factor(df_to_plot$Var1)
df_to_plot=df_to_plot[df_to_plot$Cancer %in% c('yes', 'no'),]
df_to_plot['key_names']<-df_to_plot[x_group]

df_to_plot<-relabel_objectives(df_to_plot)

#df_to_plot %>% 
 # group_by(c(Cancer, objective, key_names, Var1)) 

show_p<-plotbyObjective(df_to_plot, 'Methods' )

# TODO REGROUP THE TOTALS 

show_p

new_concise<-new_concise[!is.na(new_concise['Data']),]
new_concise<-new[c('Data', 'objective', 'method' )]







##### Section 4: PRINT METHODS to table

width=10

x_group<-'method_orig'
colname='method'
width=7


new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))
new[colname] <-apply(new[colname], 1, function(x) trimws(tolower(x)))




#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
new<-new[! (is.na(new[x_group]) | is.na(new['Data'] )),]



df_by_group_meth<-new %>%
  group_by_at(x_group) %>%
  group_modify(~ preprocessing(.x, colname) %>%
                 get_frequencies() 
  )  
df_by_group_meth<-df_by_group_meth[order(df_by_group_meth$Var1),]
aggr_methods<-aggregate(method_orig ~., df_by_group_meth[c('Var1', 'method_orig')], toString)
aggr_methods$x<-'\\\\'

write.table(df_by_group_meth, file = "review/output/methods_freq.txt", sep='\t', row.names = FALSE, quote = FALSE)
df_by_group_meth$x<-'\\\\'
write.table(df_by_group_meth, file = "review/output/methods_freq_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)

df_by_group_meth_obj<-new %>%
  group_by_at(c('Cancer', x_group, 'objective')) %>%
  group_modify(~ preprocessing(.x, colname) %>%
                 get_frequencies() 
  )  
df_by_group_meth_obj<-df_by_group_meth_obj[
  with(df_by_group_meth_obj, order(objective, Var1, Cancer, Freq, decreasing = TRUE)),
  ]
write.table(df_by_group_meth_obj, file = "review/output/objective_methods_freq.txt", sep='\t', row.names = FALSE, quote = FALSE)
df_by_group_meth_obj$x<-'\\\\'

write.table(df_by_group_meth_obj, file = "review/output/objective_methods_freq_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)

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


#new2<-new2 %>% filter(Cancer %in% cancer_filter)


new2$Data<-tolower(trimws(new2$Data))
new2<-new2 %>% filter(Data %in% tolower(level1))

levels(as.factor(new2$Data))

new2<-new2[!is.na(new2$method),]
new2<-new2[!is.na(new2$Data),]
new2<-new2[!is.na(new2$objective),]

axis1='Data'
axis2='objective'

counts<-new2 %>% count(objective, method)

if (cancer_filter == 'yes'){n_cutoff=3} else {n_cutoff=4}
counts<-counts%>% filter(n>n_cutoff)


df<-counts
counts<-new2 %>% count(objective, method)
counts<-counts%>% filter(n>n_cutoff)




# New implementation with ggalluvial
axis1='Data'
axis2="objective"
counts <- new2 %>% 
  count( Data, objective) %>% 
  mutate(
    col = objective
  ) %>%
  ggalluvial::to_lodes_form(key = type, axes = c(axis1, axis2))

df<-counts %>% filter(n>n_cutoff)
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



stats


# Extract information 
## TODO: CREATE a new column whch is merged from objective and method
# TODO: add a \newline after each objective 
  stats_to_write<-stats[c('PMID','Title', 'Data', 'ObjeMeth')]
stats_to_write$x<-'\\\\'

stats_to_write$ObjeMeth<- gsub(',',',\\\\newline',stats_to_write$ObjeMeth)
#install.packages('gdata')

write.table(stats_to_write, file = "review/output/literature_latex_2.txt", sep = " & ", row.names = FALSE, quote = FALSE)



