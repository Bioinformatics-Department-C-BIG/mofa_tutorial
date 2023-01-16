
library(ggplot2)
library(reshape2)
#install_phantomjs(version = "2.1.1",
#                  baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/",
#                  force = FALSE)
library('flextable')
#install.packages('webshot')
library('webshot')
library(dplyr)
library(tidyr)


library('readxl')
#source('literature_review_combined.R')
sh1<-read_excel(paste0(data_int_dir,'MultiOmics Tools Evaluation_2.xlsx' ),sheet=1)
colnames(sh1)
sh1<-sh1[,c(2,3:17)]
sh1<-sh1[-c(2),]


sh1 <- sh1[!((is.na(sh1[1]=='') | (sh1[1]=='') )),]
cat_names<-colnames(sh1[,c(8:16)])
sh1<-data.frame(sh1)
rownames(sh1)<-sh1$Name

# Concatenate categories first 
colnames(sh1)[1]<-'tool'
method_cats<-sh1[,c(1,8:16)]

# remove empty rownames
rownames(method_cats)<-method_cats$tool
method_cats$tool<-NULL
method_cats[method_cats=='X']<-1
method_cats[method_cats!=1]<-0
method_cats[is.na(method_cats)]<-0



##### Order by each column one after the other 

dd<-method_cats
method_cats_ordered<-dd[
  with(dd, order(dd[,1], dd[,2], dd[,3], dd[,4], dd[,5], dd[,6], dd[,7], decreasing = TRUE )),
]
cat_names
sel_obj<-which(colnames(method_cats_ordered) %in% c('FA',  "NB",  "KB",  "PR",  "jDR", "COR" ,"REG" ))
method_cats_ordered<-method_cats_ordered[,sel_obj]
method_cats_ordered

# now melt  
#colnames(method_cats_ordered)<-cat_names
method_cats_ordered<-  add_rownames(method_cats_ordered, var='tool')

# add the categories in one row
t2<-method_cats_ordered %>%
    gather(key, value, -tool)%>%
    filter(value==1) 

t3<-t2 %>%  
  select(tool, key)

t4<-plyr::ddply(t3, .(tool), colwise(paste), collapse = ",\\newline ")

t4<-plyr::ddply(t3, .(tool), colwise(paste), collapse = ", ")

# Print tools and cats in one row!!
write.csv(t4, 'review/output/tools_categories.csv')
t4


# Now add back the other info, datasets and objectives 
sh1_other<-sh1[c('tool', 'Datasets',  'Objectives' , 'Citations')]
to_keep<-sh1$tool[is.na(sh1['Reject'])]

t4_with_info<-merge(t4, sh1_other, by='tool' )
t4_with_info<-t4_with_info[match(method_cats_t$tool,t4_with_info$tool),]
t4_with_info<-t4_with_info[t4_with_info$tool %in% to_keep,]
# TODO: concat these categories with the other al_tools 



t4_with_info<-as.data.frame(t4_with_info)
t4_with_info<-t4_with_info[!is.na(t4_with_info$tool),]

t4_with_info_lt<-t4_with_info
t4_with_info_lt$x<-'\\\\'
write.table(t4_with_info_lt, file = "review/output/method_categories_concatenated_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)






sh2<-read_excel(paste0(data_int_dir,'MultiOmics Tools Evaluation_2.xlsx' ),sheet=3)

sh2<-sh2[c(2,5:10)]
sh2 <- sh2[!((is.na(sh2[1]=='') | (sh2[1]=='') )),]
sh2<-data.frame(sh2)
rownames(sh2)<-sh2$Name




sh2$Name<-NULL
sh2_chr<-sapply(sh2, as.character)
row.names(sh2_chr)<-row.names(sh2)
sh2$Ch1...Nonlinear.interactions[startsWith('X*',sh2$Ch1...Nonlinear.interactions)]
sh2_chr[which(startsWith(sh2_chr,'X'))]<-1

sh2_chr[sh2_chr!=1]<-NA


sh2_chr %>% as.data.frame() %>%
  add_rownames() %>% 
  flextable() 

#EXPORT TO LATEX
sh2_chr_latex=sh2_chr
df<- sh2_chr_latex[,-c(5)]
df[is.na(df)]<-' '
colnames(df)<-c('Non-linear interactions', 'Unequal sizes', 'Missing Data', 'NP problem', 'Heterog. Datasets')
df[df==1]<-as.character('\\checkmark')
df<-as.data.frame(df)
df_challenges<-df

df$x<-'\\\\'
write.table(df, file = "review/output/method_challenges.txt", sep=' & ', row.names = TRUE, quote = FALSE)

#ALIGN ALL TABLES 
df_challenges<-  add_rownames(df_challenges, var='tool')


t4_with_info$tool_lc<-tolower(t4_with_info$tool)
df_challenges$tool_lc<-tolower(df_challenges$tool)

t_final<-merge(t4_with_info, df_challenges, by='tool_lc' )

t_final$tool<-t_final$tool.x
t_final<-t_final[!is.na(t_final$tool),]
t_final<-t_final[match(method_cats_t$tool,t_final$tool),]

t_final<-t_final[,!names(t_final) %in% 
     c("tool.x", "tool.y", "tool_lc")]

t_final<-t_final%>%
  select(tool, everything())


t_final<-t_final[!is.na(t_final$tool),]

t_final$x<-'\\\\'

write.table(t_final, file = "review/output/method_categories_concatenated__challenges_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)

df_challenges$tool


t_final
method_cats_ordered


#####

# make it pretty!

ft3=flextable(t_final) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))
ft3<-set_header_labels(ft3, Var1 = 'Omics pair', Freq='Frequency', concat = 'Disease' ) 
ft3













