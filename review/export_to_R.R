#install_phantomjs(version = "2.1.1",
#                  baseURL = "https://github.com/wch/webshot/releases/download/v0.3.1/",
#                  force = FALSE)
library('flextable')
#install.packages('webshot')
library('webshot')

source('review/export_to_latex.R')

#### 

tool_obj[tool_obj$method_orig=='pls',]

all_tools<-read_excel(paste0(data_int_dir,'/all tools.xlsx'))
all_tools_file<-'D:/DATADRIVE/Efi Athieniti/Documents/Google Drive/PHD 2020/Literature/Data Integration/all_tools_updated.xlsx'
all_tools<-read_excel(all_tools_file)

all_tools_new<-as.data.frame(merge(all_tools[c('method_orig')], tool_obj[c('method_orig', 'objectives_concat')], by='method_orig'))

df1 = all_tools
df2=tool_obj[c('method_orig', 'Var1',  'objectives_concat')]
df2$objectives_concat_new<-df2$objectives_concat

all_tools_new= df1 %>% full_join(df2,by="method_orig")


all_tools_new2<-sapply(all_tools_new,as.character) # convert to char, save in new df

all_tools_new2[is.na(all_tools_new2)]<-""


write.table(all_tools_new2, file = "review/output/tools_objectives_concatenated.txt", sep='\t', row.names = FALSE, quote = FALSE)

# Only output the ones that I described

all_tools_to_write<-all_tools[!is.na(all_tools$Refined) & is.na(all_tools$Rejection) ,]
all_tools_to_write<-all_tools[!is.na(all_tools$Refined),]



all_tools_to_write<-all_tools_to_write[
  with(all_tools_to_write, order( Var1,objectives_concat, method_orig,decreasing = TRUE)),
]



all_tools_to_write$x<-'\\\\'


# todo: add column to all_tools ? 
if (!('Var1' %in% colnames(all_tools_to_write))){
  print('Missing required colnames')}

write.table(all_tools_to_write[c('name display','Year', 'Var1','objectives_concat','Tool/Method' ,'x' )], file = "review/output/tools_objectives_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)

data<-all_tools_to_write[c('name display', 'Year', 'Var1','objectives_concat', 'Tool/Method', 'Citation')]
data<-data[!(data$Var1 == 'multiomics pathway analysis'),]
data<-group_methods_to_short(data,'Var1')
data$Var1
# 

library(RColorBrewer)

colnames(data)<-c("Name",'Year', 'Category', "Objectives",
                       'Type', 'Reference' )

data_to_write<-data
levels(as.factor(data_to_write$Category))
data_to_write$Category = factor(data_to_write$Category, levels=levels_to_reorder[c(11,1,3,6,2,4,5,8,9,7,10)])


data_to_write<-with(data_to_write, data_to_write[order(Category),])

new1<-data_to_write[order('Category'),]
data_to_write$x<-'\\\\'
write.table(data_to_write, file = "review/output/tools_table_latex.txt", sep=' & ', row.names = FALSE, quote = FALSE)



new1


ft1=flextable(data) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))
  #  highlight(color=scales::col_factor(palette = "Paired", domain = NULL, alpha = 0.2), j=c("Var1"))

ft1 <- set_header_labels(ft1,
                           values = list('name display' = "Name",
                                         'objectives_concat' = "Objectives",
                                         Var1 = "Category", 
                                         'Tool/Method'='Type') )
ft1
save_as_html( 'Tools'= ft1, 
              path = 'plots/tools_table.html')

save_as_image( ft1, 
               path = 'plots/tools_table.png')



data2<-obj_tool_wide[,c('objective', 'Var1', 'no', 'yes')]

ft2=flextable(data2) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))%>%
  highlight(color=scales::col_factor(palette = "Paired", domain = NULL, alpha = 0.2), j=c("no", "yes"), source=c("Var1"))


ft2<-void(ft2, j='Var1', part='body')
ft2<-merge_v(ft2, j = 'objective', target = c('objective'), part = "body", combine = FALSE)

ft2
save_as_image( ft2, 
               path = 'plots/objective_tools_table.png')



##### Third table with disease 
data3=  df_nested_filtered[,-c(1,3)]
ft3=flextable(data3) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))
 # highlight(color=scales::col_factor(palette = "Paired", domain = NULL, alpha = 0.2), j=c("no", "yes"), source=c("Var1"))
ft3<-hline(ft3, i =5 )
ft3<-set_header_labels(ft3, Var1 = 'Omics pair', Freq='Frequency', concat = 'Disease' ) 
ft3
save_as_image( ft3, 
               path = 'plots/data_disease_table.png')



  data4<-aggr_methods[,1:2]
data4<-data4 %>% filter(Var1 %in% most_common_groups)%>%
  arrange(desc(Var1))

data4<-data4[with(data4, order(Var1)), ]

ft4=flextable(data4) %>% 
  autofit(add_w = 0.1,  part = c("body", "header"))
ft4<-hline(ft4, part = 'all' )


ft4<-set_header_labels(ft4, Var1 = 'Category', method_orig='Method' ) 

ft4
save_as_image( ft4, 
               path = 'plots/method_categories.png')

