library('plyr')

library('dplyr')
library('purrr')
source('literature_review.R')
source('utils.R')

##### Produces a set of  output files by counting number of studies for different categories

source('literature_review_combined.R')
# continue from part 3 lit review- prerequisite! 
#' TODO: check the rownames given by get frequencies..
#' TODO: use dplyr instead 
#' 
#' 
##### Section 4: PRINT METHODS to table
new<-new_objective_method
width=10

x_group<-'method_orig'
colname='method'
width=7


new[x_group] <-apply(new[x_group], 1, function(x) trimws(tolower(x)))
new[colname] <-apply(new[colname], 1, function(x) trimws(tolower(x)))

new<-new[! (is.na(new[x_group]) | is.na(new['Data'] )),]



df_by_group_meth<-new %>%
  group_by_at(x_group) %>%
  group_modify(~ preprocessing(.x, colname) %>%
                 get_frequencies() 
  )  
df_by_group_meth<-df_by_group_meth[order(df_by_group_meth$Var1),]
aggr_methods<-aggregate(method_orig ~., df_by_group_meth[c('Var1', 'method_orig')], toString)
write.table(aggr_methods, file = "review/output/methods_category.txt", sep='\t', row.names = FALSE, quote = FALSE)



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



##### Table 1: For each objective list ALL tools and color by category 

new_df<-df_by_group_meth_obj
# FILTER

new_df<-new_df[new_df$objective %in% most_common_objectives,]
new_df<-new_df[new_df$Var1 %in% most_common_groups,]
new_df<-new_df[new_df$Cancer %in% c('yes' , 'no'),]

new_df$method_freq<-paste(new_df$method_orig,new_df$Freq)

obj_tool<-new_df[,c('method_freq', 'Var1', 'objective', 'Cancer')] %>% 
  group_by_at(c('objective', 'Cancer')) %>% 
  nest(methods=c(method_freq))
  
obj_tool$methods_concat<-sapply(obj_tool$methods, function(x){
  paste(unlist(x), collapse=', ')
}
)
obj_tool_1<-obj_tool[,c('objective', 'Var1', 'methods_concat', 'Cancer')]
obj_tool_wide<-spread(obj_tool_1, Cancer, methods_concat) %>%
                      arrange(objective,Var1   )
                
write.table(obj_tool[,c('objective', 'Var1', 'methods_concat')], file = "review/output/objectives_methods_concat.txt", sep='\t', row.names = FALSE, quote = FALSE)
obj_tool[obj_tool$Var1=='jdr - nonlinear',]



##### Table 2: Selected popular tools and all the objectives they are used for 
# create labels 



new_df$objective<-relabel_objectives(new_df$objective)


tool_obj<-new_df[,c('method_orig', 'Var1', 'objective')] %>% 
  group_by_at('method_orig') %>% 
  nest(objectives=c(objective)) 



paste(unlist(tool_obj$objectives[[1]] ), collapse=', ')
tool_obj[tool_obj$method_orig=='pls',]


# concatenate objectives in one column 
tool_obj$objectives_concat<-sapply(tool_obj$objectives, function(x){
  paste(unlist(x), collapse=', ')
}
)


write.table(tool_obj[,c('method_orig', 'Var1', 'objectives_concat')], file = "review/output/tools_objectives.txt", sep='\t', row.names = FALSE, quote = FALSE)


resources<-read_excel(paste0(data_int_dir, 'MultiomicsResources.xlsx'))


resources$x<-'\\\\'

write.table(resources, file = "review/output/resources.txt", sep=' & ', row.names = FALSE, quote = FALSE)



t4<-read_excel(paste0(data_int_dir, 'Multiomics Tools Evaluation_2.xlsx'), sheet =3)

t4<-t4[c(1,2,3)]
t4$x<-'\\\\'

write.table(t4, file = "review/output/t4.txt", sep=' & ', row.names = FALSE, quote = FALSE)

